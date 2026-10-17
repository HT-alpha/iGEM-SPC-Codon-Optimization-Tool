//! `cf-primer-search` — search-engine layer.
//!
//! Sliding windows on the two flanks of a selected region, four pruning
//! filters, then hairpin and dimer evaluation on the survivors. The front is
//! minimized on Tm mismatch, worst-dimer penalty, GC imbalance and product-size
//! deviation. Ranking uses `cf-metrics`; thermodynamics use `cf-primer-thermo`.
//!
//! Declared dependencies: cf-primer-thermo, cf-metrics, cf-types.
//! Contract: `planning/02-schemas/pareto_front.schema.json`.

#![forbid(unsafe_code)]

use std::cmp::Ordering;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering as AtomicOrdering};
use std::time::Instant;

use cf_metrics::{fast_nondominated_sort, Ranked};

struct PrimerCounter {
    ns: AtomicU64,
    calls: AtomicU64,
}

impl PrimerCounter {
    const fn new() -> Self {
        Self {
            ns: AtomicU64::new(0),
            calls: AtomicU64::new(0),
        }
    }

    fn record(&self, started: Instant) {
        self.ns
            .fetch_add(started.elapsed().as_nanos() as u64, AtomicOrdering::Relaxed);
        self.calls.fetch_add(1, AtomicOrdering::Relaxed);
    }
}

struct PrimerProfile {
    enabled: AtomicBool,
    gc: PrimerCounter,
    homopolymer: PrimerCounter,
    tm: PrimerCounter,
    product: PrimerCounter,
    hairpin: PrimerCounter,
    dimer: PrimerCounter,
    pair_tm: PrimerCounter,
    rank: PrimerCounter,
}

static PRIMER_PROFILE: PrimerProfile = PrimerProfile {
    enabled: AtomicBool::new(false),
    gc: PrimerCounter::new(),
    homopolymer: PrimerCounter::new(),
    tm: PrimerCounter::new(),
    product: PrimerCounter::new(),
    hairpin: PrimerCounter::new(),
    dimer: PrimerCounter::new(),
    pair_tm: PrimerCounter::new(),
    rank: PrimerCounter::new(),
};

fn primer_profile_on() -> bool {
    PRIMER_PROFILE.enabled.load(AtomicOrdering::Relaxed)
}

fn arm_primer_profile() {
    if std::env::var_os("CF_PRIMER_PROFILE").is_some() {
        PRIMER_PROFILE.enabled.store(true, AtomicOrdering::Relaxed);
    }
}

fn dump_primer_profile() {
    if !primer_profile_on() {
        return;
    }
    let row = |name: &str, counter: &PrimerCounter| {
        let calls = counter.calls.load(AtomicOrdering::Relaxed);
        let ns = counter.ns.load(AtomicOrdering::Relaxed);
        format!(
            "{{\"name\":\"{name}\",\"calls\":{calls},\"seconds\":{:.6}}}",
            ns as f64 / 1e9
        )
    };
    eprintln!(
        "CF_PRIMER_PROFILE {}",
        [
            row("gc", &PRIMER_PROFILE.gc),
            row("homopolymer", &PRIMER_PROFILE.homopolymer),
            row("tm_prefilter", &PRIMER_PROFILE.tm),
            row("product", &PRIMER_PROFILE.product),
            row("hairpin", &PRIMER_PROFILE.hairpin),
            row("dimer", &PRIMER_PROFILE.dimer),
            row("pair_tm", &PRIMER_PROFILE.pair_tm),
            row("rank", &PRIMER_PROFILE.rank),
        ]
        .join(",")
    );
}
use cf_primer_thermo::{
    dg_zone, dimer, gc_clamp, gc_percent, hairpin, longest_homopolymer, perfect_duplex, DgZone,
    Parameters,
};
use serde_json::{json, Value};

#[derive(Debug)]
pub struct SearchError(pub String);

impl std::fmt::Display for SearchError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

#[derive(Debug, Clone)]
pub struct Config {
    pub tm_target: f64,
    pub product_min: usize,
    pub product_max: usize,
    pub flank: usize,
    pub primer_min: usize,
    pub primer_max: usize,
}

impl Config {
    pub fn from_json(value: &Value) -> Result<Self, SearchError> {
        let config = value.get("config").unwrap_or(value);
        let range = config
            .get("product_size_range")
            .ok_or_else(|| SearchError("CF-PRIM-006 product_size_range is required".into()))?;
        let lengths = config.get("primer_length_range");
        let flank = config
            .get("flank_length")
            .and_then(Value::as_u64)
            .unwrap_or(200) as usize;
        if flank > 300 {
            return Err(SearchError(
                "CF-PRIM-006 flank_length cannot exceed 300".into(),
            ));
        }
        Ok(Self {
            tm_target: config
                .get("tm_target")
                .and_then(Value::as_f64)
                .unwrap_or(60.0),
            product_min: range.get("min").and_then(Value::as_u64).ok_or_else(|| {
                SearchError("CF-PRIM-006 product_size_range.min is required".into())
            })? as usize,
            product_max: range.get("max").and_then(Value::as_u64).ok_or_else(|| {
                SearchError("CF-PRIM-006 product_size_range.max is required".into())
            })? as usize,
            flank,
            primer_min: lengths
                .and_then(|item| item.get("min"))
                .and_then(Value::as_u64)
                .unwrap_or(18) as usize,
            primer_max: lengths
                .and_then(|item| item.get("max"))
                .and_then(Value::as_u64)
                .unwrap_or(25) as usize,
        })
    }
}

#[derive(Debug, Clone)]
struct Primer {
    sequence: String,
    start: usize,
    end: usize,
    tm: f64,
    gc: f64,
    /// Filled after the Tm sort truncates to 40. Absent until then: hairpin
    /// and self-dimer depend only on this oligo, so they are computed once
    /// per survivor rather than once per pair.
    hairpin_dg: f64,
    hairpin_diagram: String,
    self_dimer_dg: f64,
}

pub fn gc_keeps(gc: f64) -> bool {
    (30.0..=70.0).contains(&gc)
}

pub fn homopolymer_keeps(seq: &str) -> bool {
    longest_homopolymer(seq) < 5
}

pub fn tm_prefilter_keeps(tm: f64, target: f64) -> bool {
    (tm - target).abs() <= 8.0
}

pub fn product_keeps(size: usize, min: usize, max: usize) -> bool {
    (min..=max).contains(&size)
}

fn prim001(dropped: &std::collections::BTreeMap<&str, usize>) -> String {
    let filters = [
        (
            "gc",
            "GC 30-70%",
            "Widen the flank so more base compositions are available.",
        ),
        (
            "homopolymer",
            "homopolymer (runs of 5 or longer refused)",
            "Widen the flank. This template is a long run of one base.",
        ),
        (
            "tm",
            "Tm within 8 C of the target",
            "Move the Tm target, or widen the flank.",
        ),
        ("product", "product size", "Widen the product-size range."),
    ];
    let mut best: Option<(&str, &str, &str, usize)> = None;
    for (key, label, suggestion) in filters {
        let count = dropped.get(key).copied().unwrap_or(0);
        let replace = match best {
            Some((_, _, _, previous)) => count > previous,
            None => true,
        };
        if replace {
            best = Some((key, label, suggestion, count));
        }
    }
    match best {
        Some((_, _, _, 0)) | None => {
            "CF-PRIM-001 no primer windows were generated. Widen the flank or shorten the primer."
                .into()
        }
        Some((_, label, suggestion, count)) => format!(
            "CF-PRIM-001 the {label} filter eliminated the most candidates ({count}). {suggestion}"
        ),
    }
}

fn windows(
    template: &str,
    from: usize,
    to: usize,
    primer_min: usize,
    primer_max: usize,
) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    if to <= from {
        return out;
    }
    let span = &template[from..to];
    for length in primer_min..=primer_max {
        if span.len() < length {
            continue;
        }
        for offset in 0..=span.len() - length {
            out.push((from + offset, from + offset + length));
        }
    }
    out
}

fn evaluate_primer(
    parameters: &Parameters,
    template: &str,
    start: usize,
    end: usize,
    reverse: bool,
    target: f64,
) -> Result<Result<Primer, &'static str>, SearchError> {
    let slice = &template[start..end];
    let sequence = if reverse {
        cf_primer_thermo::duplex::reverse_complement(slice).map_err(|err| SearchError(err.0))?
    } else {
        slice.to_string()
    };
    let gc_started = primer_profile_on().then(Instant::now);
    let gc_ok = gc_keeps(gc_percent(&sequence));
    if let Some(gc_started) = gc_started {
        PRIMER_PROFILE.gc.record(gc_started);
    }
    if !gc_ok {
        return Ok(Err("gc"));
    }
    let homo_started = primer_profile_on().then(Instant::now);
    let homo_ok = homopolymer_keeps(&sequence);
    if let Some(homo_started) = homo_started {
        PRIMER_PROFILE.homopolymer.record(homo_started);
    }
    if !homo_ok {
        return Ok(Err("homopolymer"));
    }
    let tm_started = primer_profile_on().then(Instant::now);
    let duplex = perfect_duplex(parameters, &sequence).map_err(|err| SearchError(err.0))?;
    if let Some(tm_started) = tm_started {
        PRIMER_PROFILE.tm.record(tm_started);
    }
    if !tm_prefilter_keeps(duplex.tm_celsius, target) {
        return Ok(Err("tm"));
    }
    Ok(Ok(Primer {
        sequence,
        start,
        end,
        tm: duplex.tm_celsius,
        gc: gc_percent(&template[start..end]),
        hairpin_dg: 0.0,
        hairpin_diagram: String::new(),
        self_dimer_dg: 0.0,
    }))
}

fn annotate_primer(parameters: &Parameters, primer: &mut Primer) -> Result<(), SearchError> {
    let hairpin_started = primer_profile_on().then(Instant::now);
    let folded = hairpin(parameters, &primer.sequence).map_err(|err| SearchError(err.0))?;
    if let Some(hairpin_started) = hairpin_started {
        PRIMER_PROFILE.hairpin.record(hairpin_started);
    }
    let dimer_started = primer_profile_on().then(Instant::now);
    let self_dimer =
        dimer(parameters, &primer.sequence, &primer.sequence).map_err(|err| SearchError(err.0))?;
    if let Some(dimer_started) = dimer_started {
        PRIMER_PROFILE.dimer.record(dimer_started);
    }
    primer.hairpin_dg = folded.delta_g_kcal;
    primer.hairpin_diagram = folded.diagram;
    primer.self_dimer_dg = self_dimer.delta_g_kcal;
    Ok(())
}

#[derive(Debug, Clone)]
struct Scored {
    objectives: [f64; 4],
    body: Value,
    refused: bool,
}

fn flag(property: &str, threshold: &str, observed: f64, zone: DgZone) -> Value {
    json!({
        "property": property,
        "threshold": threshold,
        "observed": observed,
        "severity": match zone { DgZone::Refuse => "red", DgZone::Flag => "amber", DgZone::Accept => "green" },
    })
}

fn score_pair(
    parameters: &Parameters,
    forward: &Primer,
    reverse: &Primer,
    config: &Config,
) -> Result<Option<Scored>, SearchError> {
    let product_started = primer_profile_on().then(Instant::now);
    let product = reverse.end.saturating_sub(forward.start);
    let product_ok = product_keeps(product, config.product_min, config.product_max);
    if let Some(product_started) = product_started {
        PRIMER_PROFILE.product.record(product_started);
    }
    if !product_ok {
        return Ok(None);
    }
    // Hairpin, self-dimer and Tm were measured on the oligo. Recomputing them
    // for every partner repeats the same floats. Hetero-dimer is the term
    // that depends on the pair.
    let hairpin_f_dg = forward.hairpin_dg;
    let hairpin_r_dg = reverse.hairpin_dg;
    let hairpin_f_diagram = forward.hairpin_diagram.clone();
    let hairpin_r_diagram = reverse.hairpin_diagram.clone();
    let self_f_dg = forward.self_dimer_dg;
    let self_r_dg = reverse.self_dimer_dg;
    let dimer_started = primer_profile_on().then(Instant::now);
    let hetero = dimer(parameters, &forward.sequence, &reverse.sequence)
        .map_err(|err| SearchError(err.0))?;
    if let Some(dimer_started) = dimer_started {
        PRIMER_PROFILE.dimer.record(dimer_started);
    }
    let tm_f = forward.tm;
    let tm_r = reverse.tm;
    let mismatch = (tm_f - tm_r).abs();
    let mut flags = Vec::new();
    let mut refused = false;
    for (name, value) in [
        ("hairpin_dG_forward", hairpin_f_dg),
        ("hairpin_dG_reverse", hairpin_r_dg),
        ("self_dimer_dG_forward", self_f_dg),
        ("self_dimer_dG_reverse", self_r_dg),
        ("hetero_dimer_dG", hetero.delta_g_kcal),
    ] {
        let zone = dg_zone(value);
        if zone != DgZone::Accept {
            flags.push(flag(
                name,
                if zone == DgZone::Refuse {
                    "ΔG < -9 kcal/mol"
                } else {
                    "ΔG < -4 kcal/mol"
                },
                value,
                zone,
            ));
        }
        refused |= zone == DgZone::Refuse;
    }
    let mismatch_zone = if mismatch > 5.0 {
        DgZone::Refuse
    } else if mismatch > 2.0 {
        DgZone::Flag
    } else {
        DgZone::Accept
    };
    if mismatch_zone != DgZone::Accept {
        flags.push(flag(
            "tm_mismatch",
            if mismatch_zone == DgZone::Refuse {
                "Tm mismatch > 5 C"
            } else {
                "Tm mismatch > 2 C"
            },
            mismatch,
            mismatch_zone,
        ));
    }
    refused |= mismatch_zone == DgZone::Refuse;
    let worst = [
        hairpin_f_dg,
        hairpin_r_dg,
        self_f_dg,
        self_r_dg,
        hetero.delta_g_kcal,
    ]
    .into_iter()
    .fold(0.0, f64::min);
    let midpoint = (config.product_min + config.product_max) as f64 / 2.0;
    let objectives = [
        mismatch,
        -worst,
        (forward.gc - reverse.gc).abs(),
        (product as f64 - midpoint).abs(),
    ];
    // The amplicon string is attached to the pairs that are emitted. Building
    // it here would allocate one for every survivor, and the rank keeps 20.
    let body = json!({
        "forward_sequence": forward.sequence,
        "reverse_sequence": reverse.sequence,
        "forward_start": forward.start,
        "forward_end": forward.end,
        "reverse_start": reverse.start,
        "reverse_end": reverse.end,
        "product_size": product,
        "tm_forward": tm_f,
        "tm_reverse": tm_r,
        "tm_mismatch": mismatch,
        "hairpin_dG_forward": hairpin_f_dg,
        "hairpin_dG_reverse": hairpin_r_dg,
        "hairpin_diagram_forward": hairpin_f_diagram,
        "hairpin_diagram_reverse": hairpin_r_diagram,
        "self_dimer_dG_forward": self_f_dg,
        "self_dimer_dG_reverse": self_r_dg,
        "hetero_dimer_dG": hetero.delta_g_kcal,
        "hetero_dimer_3prime_dG": hetero.three_prime_delta_g_kcal,
        "hetero_dimer_alignment": hetero.alignment,
        "gc_forward": forward.gc,
        "gc_reverse": reverse.gc,
        "gc_clamp_forward": gc_clamp(&forward.sequence),
        "gc_clamp_reverse": gc_clamp(&reverse.sequence),
        "flags": flags,
        "refused": refused,
        "parameter_set_id": parameters.parameter_set_id,
        "parameter_digest": parameters.digest,
        "parameter_citation": parameters.citation,
    });
    Ok(Some(Scored {
        objectives,
        body,
        refused,
    }))
}

pub fn analyze_pair(forward: &str, reverse: &str) -> Result<Value, SearchError> {
    let parameters = cf_primer_thermo::bundled().map_err(|err| SearchError(err.0))?;
    // A manual pair has no amplicon. Product size stays null rather than a made-up length.
    let hairpin_f = hairpin(&parameters, forward).map_err(|err| SearchError(err.0))?;
    let hairpin_r = hairpin(&parameters, reverse).map_err(|err| SearchError(err.0))?;
    let self_f = dimer(&parameters, forward, forward).map_err(|err| SearchError(err.0))?;
    let self_r = dimer(&parameters, reverse, reverse).map_err(|err| SearchError(err.0))?;
    let hetero = dimer(&parameters, forward, reverse).map_err(|err| SearchError(err.0))?;
    let tm_f = perfect_duplex(&parameters, forward).map_err(|err| SearchError(err.0))?;
    let tm_r = perfect_duplex(&parameters, reverse).map_err(|err| SearchError(err.0))?;
    let mismatch = (tm_f.tm_celsius - tm_r.tm_celsius).abs();
    let mut flags = Vec::new();
    let mut refused = false;
    for (name, value) in [
        ("hairpin_dG_forward", hairpin_f.delta_g_kcal),
        ("hairpin_dG_reverse", hairpin_r.delta_g_kcal),
        ("self_dimer_dG_forward", self_f.delta_g_kcal),
        ("self_dimer_dG_reverse", self_r.delta_g_kcal),
        ("hetero_dimer_dG", hetero.delta_g_kcal),
    ] {
        let zone = dg_zone(value);
        if zone != DgZone::Accept {
            flags.push(flag(
                name,
                if zone == DgZone::Refuse {
                    "ΔG < -9 kcal/mol"
                } else {
                    "ΔG < -4 kcal/mol"
                },
                value,
                zone,
            ));
        }
        refused |= zone == DgZone::Refuse;
    }
    if mismatch > 5.0 {
        flags.push(flag(
            "tm_mismatch",
            "Tm mismatch > 5 C",
            mismatch,
            DgZone::Refuse,
        ));
        refused = true;
    } else if mismatch > 2.0 {
        flags.push(flag(
            "tm_mismatch",
            "Tm mismatch > 2 C",
            mismatch,
            DgZone::Flag,
        ));
    }
    Ok(json!({
        "forward_sequence": forward,
        "reverse_sequence": reverse,
        "product_size": null,
        "product_sequence": null,
        "tm_forward": tm_f.tm_celsius,
        "tm_reverse": tm_r.tm_celsius,
        "tm_mismatch": mismatch,
        "hairpin_dG_forward": hairpin_f.delta_g_kcal,
        "hairpin_dG_reverse": hairpin_r.delta_g_kcal,
        "hairpin_diagram_forward": hairpin_f.diagram,
        "hairpin_diagram_reverse": hairpin_r.diagram,
        "self_dimer_dG_forward": self_f.delta_g_kcal,
        "self_dimer_dG_reverse": self_r.delta_g_kcal,
        "hetero_dimer_dG": hetero.delta_g_kcal,
        "hetero_dimer_3prime_dG": hetero.three_prime_delta_g_kcal,
        "hetero_dimer_alignment": hetero.alignment,
        "gc_forward": gc_percent(forward),
        "gc_reverse": gc_percent(reverse),
        "gc_clamp_forward": gc_clamp(forward),
        "gc_clamp_reverse": gc_clamp(reverse),
        "flags": flags,
        "refused": refused,
        "parameter_set_id": parameters.parameter_set_id,
        "parameter_digest": parameters.digest,
        "parameter_citation": parameters.citation,
        "blast": "Specificity against genome not checked.",
    }))
}

pub fn search(
    template: &str,
    start: usize,
    end: usize,
    config: &Config,
    progress: &mut dyn FnMut(f64),
) -> Result<Value, SearchError> {
    arm_primer_profile();
    let parameters = cf_primer_thermo::bundled().map_err(|err| SearchError(err.0))?;
    progress(0.02);
    if end > template.len() || start >= end {
        return Err(SearchError(format!(
            "CF-PRIM-003 region {start}..{end} is outside the sequence of length {}",
            template.len()
        )));
    }
    let upstream_from = start.saturating_sub(config.flank);
    let downstream_to = (end + config.flank).min(template.len());
    if start - upstream_from < config.primer_min || downstream_to - end < config.primer_min {
        return Err(SearchError(
            "CF-PRIM-002 the sequence is too short for a primer on the requested flank".into(),
        ));
    }
    let mut dropped = std::collections::BTreeMap::<&str, usize>::new();
    let mut note_drop = |reason: &'static str| {
        *dropped.entry(reason).or_insert(0) += 1;
    };
    let mut forward = Vec::new();
    for (window_start, window_end) in windows(
        template,
        upstream_from,
        start,
        config.primer_min,
        config.primer_max,
    ) {
        match evaluate_primer(
            &parameters,
            template,
            window_start,
            window_end,
            false,
            config.tm_target,
        )? {
            Ok(primer) => forward.push(primer),
            Err(reason) => note_drop(reason),
        }
    }
    let mut reverse = Vec::new();
    for (window_start, window_end) in windows(
        template,
        end,
        downstream_to,
        config.primer_min,
        config.primer_max,
    ) {
        match evaluate_primer(
            &parameters,
            template,
            window_start,
            window_end,
            true,
            config.tm_target,
        )? {
            Ok(primer) => reverse.push(primer),
            Err(reason) => note_drop(reason),
        }
    }
    progress(0.25);
    forward.sort_by(|left, right| {
        (left.tm - config.tm_target)
            .abs()
            .partial_cmp(&(right.tm - config.tm_target).abs())
            .unwrap_or(Ordering::Equal)
    });
    reverse.sort_by(|left, right| {
        (left.tm - config.tm_target)
            .abs()
            .partial_cmp(&(right.tm - config.tm_target).abs())
            .unwrap_or(Ordering::Equal)
    });
    let mut truncated = false;
    if forward.len() > 40 {
        forward.truncate(40);
        truncated = true;
    }
    if reverse.len() > 40 {
        reverse.truncate(40);
        truncated = true;
    }
    for primer in forward.iter_mut().chain(reverse.iter_mut()) {
        annotate_primer(&parameters, primer)?;
    }
    let total = (forward.len() * reverse.len()).max(1);
    let mut scored = Vec::new();
    for (index, fwd) in forward.iter().enumerate() {
        for rev in &reverse {
            match score_pair(&parameters, fwd, rev, config)? {
                Some(pair) => scored.push(pair),
                None => note_drop("product"),
            }
        }
        progress(0.25 + 0.7 * ((index + 1) * reverse.len()) as f64 / total as f64);
    }
    if scored.is_empty() {
        dump_primer_profile();
        return Err(SearchError(prim001(&dropped)));
    }
    let rank_started = primer_profile_on().then(Instant::now);
    let ranked: Vec<Ranked> = scored
        .iter()
        .map(|pair| {
            let objectives = pair.objectives.iter().copied().map(Some).collect();
            if pair.refused {
                Ranked::infeasible(objectives, 1.0)
            } else {
                Ranked::feasible(objectives)
            }
        })
        .collect();
    let fronts = fast_nondominated_sort(&ranked);
    let front = fronts.first().cloned().unwrap_or_default();
    let mut ordered = front;
    ordered.sort_by(|left, right| {
        scored[*left].objectives[0]
            .partial_cmp(&scored[*right].objectives[0])
            .unwrap_or(Ordering::Equal)
            .then_with(|| {
                scored[*left].objectives[1]
                    .partial_cmp(&scored[*right].objectives[1])
                    .unwrap_or(Ordering::Equal)
            })
    });
    ordered.truncate(20);
    let recommended = ordered.first().copied();
    let pairs: Vec<Value> = ordered
        .iter()
        .enumerate()
        .map(|(rank_index, index)| {
            let mut body = scored[*index].body.clone();
            let start = body["forward_start"].as_u64().expect("forward_start") as usize;
            let end = body["reverse_end"].as_u64().expect("reverse_end") as usize;
            body["product_sequence"] = json!(&template[start..end]);
            body["pair_id"] = json!(format!("pair_{rank_index:02}_{:08x}", *index));
            body["nondomination_rank"] = json!(0);
            body["objective_vector"] = json!(scored[*index].objectives);
            body
        })
        .collect();
    if let Some(rank_started) = rank_started {
        PRIMER_PROFILE.rank.record(rank_started);
    }
    progress(1.0);
    dump_primer_profile();
    Ok(json!({
        "parameter_set_id": parameters.parameter_set_id,
        "parameter_digest": parameters.digest,
        "parameter_citation": parameters.citation,
        "objective_order": [
            "primer.tm_mismatch",
            "primer.worst_dimer_penalty",
            "primer.gc_balance",
            "primer.product_deviation"
        ],
        "objective_sense": ["minimize", "minimize", "minimize", "minimize"],
        "pairs": pairs,
        "recommended_index": recommended.map(|_| 0),
        "survivors": scored.len(),
        "truncated": truncated,
        "blast": "Specificity against genome not checked.",
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn each_prune_filter_removes_what_it_should() {
        assert!(!gc_keeps(80.0));
        assert!(gc_keeps(50.0));
        assert!(!homopolymer_keeps("ACGAAAAAT"));
        assert!(homopolymer_keeps("ACGTACGT"));
        assert!(!tm_prefilter_keeps(40.0, 60.0));
        assert!(tm_prefilter_keeps(55.0, 60.0));
        assert!(!product_keeps(50, 100, 200));
        assert!(product_keeps(150, 100, 200));
    }

    #[test]
    fn a_short_flank_is_a_named_refusal() {
        let config = Config {
            tm_target: 60.0,
            product_min: 40,
            product_max: 80,
            flank: 200,
            primer_min: 18,
            primer_max: 25,
        };
        let error = search("ACGTACGT", 4, 6, &config, &mut |_| {}).unwrap_err();
        assert!(error.0.contains("CF-PRIM-002"), "{error}");
    }

    #[test]
    fn a_homopolymer_template_names_that_filter() {
        let config = Config {
            tm_target: 60.0,
            product_min: 20,
            product_max: 80,
            flank: 40,
            primer_min: 18,
            primer_max: 20,
        };
        let template = "AAACCGGTAA".repeat(20);
        let error = search(&template, 40, 60, &config, &mut |_| {}).unwrap_err();
        assert!(error.0.contains("homopolymer"), "{error}");
        assert!(error.0.contains("CF-PRIM-001"), "{error}");
    }

    #[test]
    fn emitted_pairs_embed_the_template_slice() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../benchmarks/validation/sequences/rpsL.fa"
        );
        let template: String = std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('>'))
            .collect();
        let flank = 40;
        let start = flank;
        let end = template.len() - flank;
        let region = end - start;
        let config = Config {
            tm_target: 60.0,
            product_min: region,
            product_max: region + 2 * flank,
            flank,
            primer_min: 18,
            primer_max: 25,
        };
        let document = search(&template, start, end, &config, &mut |_| {}).unwrap();
        let pairs = document["pairs"].as_array().unwrap();
        assert!(!pairs.is_empty());
        let survivors = document["survivors"].as_u64().unwrap();
        assert!(survivors >= pairs.len() as u64);
        for pair in pairs {
            let from = pair["forward_start"].as_u64().unwrap() as usize;
            let to = pair["reverse_end"].as_u64().unwrap() as usize;
            let product = pair["product_sequence"].as_str().unwrap();
            assert_eq!(product, &template[from..to]);
            assert_eq!(
                product.len(),
                pair["product_size"].as_u64().unwrap() as usize
            );
        }
    }

    #[test]
    fn search_returns_a_front_or_cf_prim_001() {
        let template = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG";
        let config = Config {
            tm_target: 60.0,
            product_min: 20,
            product_max: 90,
            flank: 40,
            primer_min: 18,
            primer_max: 22,
        };
        match search(template, 30, 55, &config, &mut |_| {}) {
            Ok(document) => {
                assert!(!document["pairs"].as_array().unwrap().is_empty());
                assert_eq!(document["objective_order"].as_array().unwrap().len(), 4);
            }
            Err(error) => assert!(error.0.contains("CF-PRIM-001"), "{error}"),
        }
    }
}
