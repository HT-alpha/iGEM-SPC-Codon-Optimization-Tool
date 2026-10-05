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

use cf_metrics::{fast_nondominated_sort, Ranked};
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
) -> Result<Option<Primer>, SearchError> {
    let slice = &template[start..end];
    let sequence = if reverse {
        cf_primer_thermo::duplex::reverse_complement(slice).map_err(|err| SearchError(err.0))?
    } else {
        slice.to_string()
    };
    if !gc_keeps(gc_percent(&sequence)) || !homopolymer_keeps(&sequence) {
        return Ok(None);
    }
    let duplex = perfect_duplex(parameters, &sequence).map_err(|err| SearchError(err.0))?;
    if !tm_prefilter_keeps(duplex.tm_celsius, target) {
        return Ok(None);
    }
    Ok(Some(Primer {
        sequence,
        start,
        end,
        tm: duplex.tm_celsius,
        gc: gc_percent(&template[start..end]),
    }))
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
    template: &str,
    forward: &Primer,
    reverse: &Primer,
    config: &Config,
) -> Result<Option<Scored>, SearchError> {
    let product = reverse.end.saturating_sub(forward.start);
    if !product_keeps(product, config.product_min, config.product_max) {
        return Ok(None);
    }
    let hairpin_f = hairpin(parameters, &forward.sequence).map_err(|err| SearchError(err.0))?;
    let hairpin_r = hairpin(parameters, &reverse.sequence).map_err(|err| SearchError(err.0))?;
    let self_f = dimer(parameters, &forward.sequence, &forward.sequence)
        .map_err(|err| SearchError(err.0))?;
    let self_r = dimer(parameters, &reverse.sequence, &reverse.sequence)
        .map_err(|err| SearchError(err.0))?;
    let hetero = dimer(parameters, &forward.sequence, &reverse.sequence)
        .map_err(|err| SearchError(err.0))?;
    let tm_f = perfect_duplex(parameters, &forward.sequence).map_err(|err| SearchError(err.0))?;
    let tm_r = perfect_duplex(parameters, &reverse.sequence).map_err(|err| SearchError(err.0))?;
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
        hairpin_f.delta_g_kcal,
        hairpin_r.delta_g_kcal,
        self_f.delta_g_kcal,
        self_r.delta_g_kcal,
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
    let product_sequence = template[forward.start..reverse.end].to_string();
    let body = json!({
        "forward_sequence": forward.sequence,
        "reverse_sequence": reverse.sequence,
        "forward_start": forward.start,
        "forward_end": forward.end,
        "reverse_start": reverse.start,
        "reverse_end": reverse.end,
        "product_size": product,
        "product_sequence": product_sequence,
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
    let mut forward = Vec::new();
    for (window_start, window_end) in windows(
        template,
        upstream_from,
        start,
        config.primer_min,
        config.primer_max,
    ) {
        if let Some(primer) = evaluate_primer(
            &parameters,
            template,
            window_start,
            window_end,
            false,
            config.tm_target,
        )? {
            forward.push(primer);
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
        if let Some(primer) = evaluate_primer(
            &parameters,
            template,
            window_start,
            window_end,
            true,
            config.tm_target,
        )? {
            reverse.push(primer);
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
    let total = (forward.len() * reverse.len()).max(1);
    let mut scored = Vec::new();
    for (index, fwd) in forward.iter().enumerate() {
        for rev in &reverse {
            if let Some(pair) = score_pair(&parameters, template, fwd, rev, config)? {
                scored.push(pair);
            }
        }
        progress(0.25 + 0.7 * ((index + 1) * reverse.len()) as f64 / total as f64);
    }
    if scored.is_empty() {
        return Err(SearchError(
            "CF-PRIM-001 no primer pairs survived GC, homopolymer, Tm and product-size pruning. Widen the flank, widen the product-size range, or move the Tm target.".into(),
        ));
    }
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
            body["pair_id"] = json!(format!("pair_{rank_index:02}_{:08x}", *index));
            body["nondomination_rank"] = json!(0);
            body["objective_vector"] = json!(scored[*index].objectives);
            body
        })
        .collect();
    progress(1.0);
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
