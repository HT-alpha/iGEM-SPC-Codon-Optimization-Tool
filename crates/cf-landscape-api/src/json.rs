//! Document form for the query and response types.
//!
//! Ticket P1-20, `TST-LAPI-001`.
//!
//! ADR-0001 accepts type duplication across the Rust/Python boundary as the
//! price of a Rust hot path. The bill for that acceptance is this module plus
//! the contract test that drives it: a document written by Rust, re-read by
//! Rust, handed to Python and canonicalized there must come back byte for byte.
//! Drift between the two views of the contract then shows up as a test failure
//! rather than as a field that silently means something different on each side.
//!
//! ## Two deliberate asymmetries
//!
//! **A `gradients` block is refused on parse.** The response schema admits one
//! because a provider may offer gradients to some consumer outside DNA Design Suite.
//! No query this crate can construct asks for one — [`Outputs`] has no such
//! field — so a response carrying gradients is answering a question this
//! contract did not pose. Parsing it would mean inventing a semantics for a
//! path ADR-0002 closed; refusing names the decision instead.
//!
//! **`construct_context` round-trips without being interpreted.** Its inline
//! form is a large schema owned by the ingest layer, and no Phase 1 term reads
//! it. Preserving it verbatim is honest about that: the document survives, and
//! nothing here pretends to understand it.

use cf_types::{
    canonicalize, code, ArtifactId, CfError, Codon, CodonSequence, Digest, Direction, HostId,
    Result, TermId, Unit,
};
use serde_json::{json, Map, Value};

use crate::capability::Capability;
use crate::query::{
    Budget, BudgetTier, CachePolicy, Edit, HostContext, LandscapeQuery, ModelPin, Outputs, Target,
};
use crate::response::{
    AlternativeDelta, DeltaEntry, Deltas, Diagnostics, LandscapeResponse, ObjectiveVector,
    Provenance, ResolvedBy, Status, TermValue, Uncertainty, UncertaintyMethod,
};
use crate::role::{ObjectiveAxis, TermRole};
use crate::{CONTRACT_MAJOR, SCHEMA_VERSION};

fn bad(detail: impl Into<String>) -> CfError {
    CfError::new(code::LAND_DOCUMENT_INVALID, detail.into())
}

fn obj<'a>(value: &'a Value, field: &str) -> Result<&'a Map<String, Value>> {
    value
        .as_object()
        .ok_or_else(|| bad(format!("{field} must be an object")))
}

fn get<'a>(map: &'a Map<String, Value>, field: &str) -> Result<&'a Value> {
    map.get(field)
        .ok_or_else(|| bad(format!("missing required field {field}")))
}

fn as_str<'a>(map: &'a Map<String, Value>, field: &str) -> Result<&'a str> {
    get(map, field)?
        .as_str()
        .ok_or_else(|| bad(format!("{field} must be a string")))
}

fn as_u64(map: &Map<String, Value>, field: &str) -> Result<u64> {
    get(map, field)?
        .as_u64()
        .ok_or_else(|| bad(format!("{field} must be a non-negative integer")))
}

fn as_f64(map: &Map<String, Value>, field: &str) -> Result<f64> {
    get(map, field)?
        .as_f64()
        .ok_or_else(|| bad(format!("{field} must be a number")))
}

fn opt_str<'a>(map: &'a Map<String, Value>, field: &str) -> Option<&'a str> {
    map.get(field).and_then(Value::as_str)
}

fn opt_bool(map: &Map<String, Value>, field: &str, default: bool) -> bool {
    map.get(field).and_then(Value::as_bool).unwrap_or(default)
}

fn digest(map: &Map<String, Value>, field: &str) -> Result<Digest> {
    Digest::parse(as_str(map, field)?)
}

fn envelope(map: &Map<String, Value>, kind: &str) -> Result<()> {
    let major = as_u64(map, "contract_major")?;
    if major != u64::from(CONTRACT_MAJOR) {
        return Err(CfError::new(
            code::RESOLVE_MAJOR_MISMATCH,
            format!("{kind} declares contract_major {major}; this build speaks {CONTRACT_MAJOR}"),
        ));
    }
    let version = as_str(map, "schema_version")?;
    if version.split('.').next() != SCHEMA_VERSION.split('.').next() {
        return Err(bad(format!(
            "{kind} declares schema_version {version}, incompatible with {SCHEMA_VERSION}"
        )));
    }
    Ok(())
}

/// A construct context carried through the contract without interpretation.
#[derive(Clone, PartialEq, Debug)]
pub enum ConstructRef {
    /// `{"kind": "ref", "construct_digest": "sha256:..."}`.
    Digest(Digest),
    /// The inline document, preserved verbatim.
    Inline(Value),
}

// --------------------------------------------------------------------------- //
// query
// --------------------------------------------------------------------------- //

impl Target {
    fn to_json(&self) -> Result<Value> {
        let mut map = Map::new();
        map.insert("kind".into(), json!(self.kind()));
        match self {
            Target::Full { sequence, .. } => {
                map.insert("sequence".into(), sequence.to_json()?);
            }
            Target::Delta {
                base_digest, edits, ..
            } => {
                map.insert("base_digest".into(), json!(base_digest.as_str()));
                map.insert(
                    "edits".into(),
                    Value::Array(
                        edits
                            .iter()
                            .map(|e| {
                                json!({
                                    "codon_index": e.codon_index,
                                    "codon": e.codon.to_string_upper(),
                                })
                            })
                            .collect(),
                    ),
                );
            }
            Target::DeltaTable {
                base_digest,
                positions,
                ..
            } => {
                map.insert("base_digest".into(), json!(base_digest.as_str()));
                map.insert("positions".into(), json!(positions));
            }
        }
        if let Some(cistron) = self.cistron_id() {
            map.insert("cistron_id".into(), json!(cistron));
        }
        Ok(Value::Object(map))
    }

    fn from_json(value: &Value) -> Result<Self> {
        let map = obj(value, "target")?;
        let cistron_id = opt_str(map, "cistron_id").map(str::to_owned);
        match as_str(map, "kind")? {
            "full" => Ok(Target::Full {
                sequence: Box::new(CodonSequence::from_json(get(map, "sequence")?)?),
                cistron_id,
            }),
            "delta" => {
                let edits = get(map, "edits")?
                    .as_array()
                    .ok_or_else(|| bad("target.edits must be an array"))?
                    .iter()
                    .map(|e| {
                        let e = obj(e, "target.edits[]")?;
                        Ok(Edit::new(
                            as_u64(e, "codon_index")? as usize,
                            Codon::parse(as_str(e, "codon")?)?,
                        ))
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(Target::Delta {
                    base_digest: digest(map, "base_digest")?,
                    edits,
                    cistron_id,
                })
            }
            "delta_table" => {
                let positions = get(map, "positions")?
                    .as_array()
                    .ok_or_else(|| bad("target.positions must be an array"))?
                    .iter()
                    .map(|p| {
                        p.as_u64()
                            .map(|v| v as usize)
                            .ok_or_else(|| bad("target.positions[] must be a non-negative integer"))
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(Target::DeltaTable {
                    base_digest: digest(map, "base_digest")?,
                    positions,
                    cistron_id,
                })
            }
            other => Err(bad(format!("unknown target kind {other:?}"))),
        }
    }
}

impl LandscapeQuery {
    /// Serialize to a `landscape_query.schema.json` document.
    pub fn to_json(&self) -> Result<Value> {
        let mut map = Map::new();
        map.insert("schema_version".into(), json!(SCHEMA_VERSION));
        map.insert("contract_major".into(), json!(CONTRACT_MAJOR));
        map.insert("query_id".into(), json!(self.query_id));
        map.insert("target".into(), self.target.to_json()?);

        let mut host = Map::new();
        host.insert("host_id".into(), json!(self.host_context.host_id.as_str()));
        host.insert(
            "host_profile_digest".into(),
            json!(self.host_context.host_profile_digest.as_str()),
        );
        if let Some(condition) = &self.host_context.condition_id {
            host.insert("condition".into(), json!({ "condition_id": condition }));
        }
        map.insert("host_context".into(), Value::Object(host));

        if let Some(construct) = &self.construct_context {
            map.insert(
                "construct_context".into(),
                match construct {
                    ConstructRef::Digest(d) => {
                        json!({ "kind": "ref", "construct_digest": d.as_str() })
                    }
                    ConstructRef::Inline(v) => v.clone(),
                },
            );
        }

        map.insert(
            "requested_terms".into(),
            Value::Array(
                self.requested_terms
                    .iter()
                    .map(|t| json!(t.as_str()))
                    .collect(),
            ),
        );
        map.insert(
            "outputs".into(),
            json!({
                "scalars": self.outputs.scalars,
                "per_position": self.outputs.per_position,
                "deltas": self.outputs.deltas,
                "uncertainty": self.outputs.uncertainty,
                "violations": self.outputs.violations,
            }),
        );

        let mut budget = Map::new();
        budget.insert("tier".into(), json!(self.budget.tier.as_str()));
        if let Some(ms) = self.budget.max_latency_ms {
            budget.insert("max_latency_ms".into(), json!(ms));
        }
        if let Some(mb) = self.budget.max_memory_mb {
            budget.insert("max_memory_mb".into(), json!(mb));
        }
        map.insert("budget".into(), Value::Object(budget));

        let mut cache = Map::new();
        cache.insert("read".into(), json!(self.cache_policy.read));
        cache.insert("write".into(), json!(self.cache_policy.write));
        if let Some(ns) = &self.cache_policy.namespace {
            cache.insert("namespace".into(), json!(ns));
        }
        map.insert("cache_policy".into(), Value::Object(cache));

        if !self.model_pins.is_empty() {
            let mut pins = Map::new();
            for (term, pin) in &self.model_pins {
                pins.insert(
                    term.as_str().to_owned(),
                    json!({
                        "artifact_id": pin.artifact_id.as_str(),
                        "digest": pin.digest.as_str(),
                    }),
                );
            }
            map.insert("model_pins".into(), Value::Object(pins));
        }

        if let Some(seed) = self.seed {
            map.insert("seed".into(), json!(seed));
        }
        Ok(Value::Object(map))
    }

    /// Parse a `landscape_query.schema.json` document, re-checking structure.
    pub fn from_json(value: &Value) -> Result<Self> {
        let map = obj(value, "query")?;
        envelope(map, "query")?;

        if map.contains_key("outputs") {
            let outputs = obj(get(map, "outputs")?, "outputs")?;
            if outputs.get("gradients").and_then(Value::as_bool) == Some(true) {
                return Err(CfError::new(
                    code::LAND_CAPABILITY,
                    "outputs.gradients is not reachable through this contract; \
                     ADR-0002 makes L2 optional and unused, and exact deltas supersede it",
                ));
            }
        }

        let host_map = obj(get(map, "host_context")?, "host_context")?;
        let mut host = HostContext::new(
            HostId::new(as_str(host_map, "host_id")?)?,
            digest(host_map, "host_profile_digest")?,
        );
        if let Some(condition) = host_map.get("condition") {
            let condition = obj(condition, "host_context.condition")?;
            host.condition_id = Some(as_str(condition, "condition_id")?.to_owned());
        }

        let requested_terms = get(map, "requested_terms")?
            .as_array()
            .ok_or_else(|| bad("requested_terms must be an array"))?
            .iter()
            .map(|t| {
                TermId::new(
                    t.as_str()
                        .ok_or_else(|| bad("requested_terms[] must be a string"))?,
                )
            })
            .collect::<Result<Vec<_>>>()?;

        let budget_map = obj(get(map, "budget")?, "budget")?;
        let mut budget = Budget::new(
            BudgetTier::parse(as_str(budget_map, "tier")?)
                .ok_or_else(|| bad("unknown budget.tier"))?,
        );
        budget.max_latency_ms = budget_map.get("max_latency_ms").and_then(Value::as_u64);
        budget.max_memory_mb = budget_map.get("max_memory_mb").and_then(Value::as_u64);

        let outputs = match map.get("outputs") {
            None => Outputs::default(),
            Some(v) => {
                let o = obj(v, "outputs")?;
                Outputs {
                    scalars: opt_bool(o, "scalars", true),
                    per_position: opt_bool(o, "per_position", false),
                    deltas: opt_bool(o, "deltas", false),
                    uncertainty: opt_bool(o, "uncertainty", false),
                    violations: opt_bool(o, "violations", false),
                }
            }
        };

        let cache_policy = match map.get("cache_policy") {
            None => CachePolicy::default(),
            Some(v) => {
                let c = obj(v, "cache_policy")?;
                CachePolicy {
                    read: opt_bool(c, "read", true),
                    write: opt_bool(c, "write", true),
                    namespace: opt_str(c, "namespace").map(str::to_owned),
                }
            }
        };

        let mut model_pins = Vec::new();
        if let Some(pins) = map.get("model_pins") {
            for (term, pin) in obj(pins, "model_pins")? {
                let pin = obj(pin, "model_pins[]")?;
                model_pins.push((
                    TermId::new(term.as_str())?,
                    ModelPin {
                        artifact_id: ArtifactId::new(as_str(pin, "artifact_id")?)?,
                        digest: digest(pin, "digest")?,
                    },
                ));
            }
            // serde_json preserves insertion order only with the preserve_order
            // feature, which is off. Sort so the parse is deterministic.
            model_pins.sort_by(|a, b| a.0.as_str().cmp(b.0.as_str()));
        }

        let construct_context = map.get("construct_context").map(|v| {
            match v
                .as_object()
                .and_then(|m| m.get("kind"))
                .and_then(Value::as_str)
            {
                Some("ref") => v
                    .as_object()
                    .and_then(|m| m.get("construct_digest"))
                    .and_then(Value::as_str)
                    .and_then(|d| Digest::parse(d).ok())
                    .map(ConstructRef::Digest)
                    .unwrap_or_else(|| ConstructRef::Inline(v.clone())),
                _ => ConstructRef::Inline(v.clone()),
            }
        });

        let query = Self {
            query_id: as_str(map, "query_id")?.to_owned(),
            target: Target::from_json(get(map, "target")?)?,
            host_context: host,
            construct_context,
            requested_terms,
            outputs,
            budget,
            cache_policy,
            model_pins,
            seed: map.get("seed").and_then(Value::as_u64),
        };
        query.validate()?;
        Ok(query)
    }

    /// RFC 8785 canonical form of the document.
    pub fn to_canonical_json(&self) -> Result<String> {
        canonicalize(&self.to_json()?)
    }
}

// --------------------------------------------------------------------------- //
// response
// --------------------------------------------------------------------------- //

fn uncertainty_to_json(u: &Uncertainty) -> Value {
    let mut map = Map::new();
    map.insert("method".into(), json!(u.method.as_str()));
    map.insert("calibrated".into(), json!(u.calibrated));
    if let Some(sigma) = u.sigma {
        map.insert("sigma".into(), json!(sigma));
    }
    if let Some((lo, hi)) = u.interval {
        map.insert("interval".into(), json!([lo, hi]));
    }
    if let Some(coverage) = u.coverage {
        map.insert("coverage".into(), json!(coverage));
    }
    if let Some(n) = u.n_members {
        map.insert("n_members".into(), json!(n));
    }
    if let Some(f) = u.shared_fraction {
        map.insert("shared_fraction".into(), json!(f));
    }
    Value::Object(map)
}

fn uncertainty_from_json(value: &Value) -> Result<Uncertainty> {
    let map = obj(value, "uncertainty")?;
    let method = UncertaintyMethod::parse(as_str(map, "method")?)
        .ok_or_else(|| bad("unknown uncertainty.method"))?;
    let interval = match map.get("interval") {
        None => None,
        Some(v) => {
            let a = v
                .as_array()
                .ok_or_else(|| bad("uncertainty.interval must be an array"))?;
            if a.len() != 2 {
                return Err(bad("uncertainty.interval must have exactly two entries"));
            }
            Some((
                a[0].as_f64().ok_or_else(|| bad("interval[0]"))?,
                a[1].as_f64().ok_or_else(|| bad("interval[1]"))?,
            ))
        }
    };
    Ok(Uncertainty {
        method,
        calibrated: get(map, "calibrated")?
            .as_bool()
            .ok_or_else(|| bad("uncertainty.calibrated must be a boolean"))?,
        sigma: map.get("sigma").and_then(Value::as_f64),
        interval,
        coverage: map.get("coverage").and_then(Value::as_f64),
        n_members: map
            .get("n_members")
            .and_then(Value::as_u64)
            .map(|n| n as u32),
        shared_fraction: map.get("shared_fraction").and_then(Value::as_f64),
    })
}

fn provenance_to_json(p: &Provenance) -> Value {
    let mut map = Map::new();
    map.insert("resolved_by".into(), json!(p.resolved_by.as_str()));
    if let Some(id) = &p.artifact_id {
        map.insert("artifact_id".into(), json!(id.as_str()));
    }
    if let Some(d) = &p.digest {
        map.insert("digest".into(), json!(d.as_str()));
    }
    if let Some(engine) = &p.engine {
        map.insert("engine".into(), json!(engine));
    }
    if let Some(version) = &p.engine_version {
        map.insert("engine_version".into(), json!(version));
    }
    if p.fallback {
        map.insert("fallback".into(), json!(true));
    }
    Value::Object(map)
}

fn provenance_from_json(value: &Value) -> Result<Provenance> {
    let map = obj(value, "provenance")?;
    Ok(Provenance {
        resolved_by: ResolvedBy::parse(as_str(map, "resolved_by")?)
            .ok_or_else(|| bad("unknown provenance.resolved_by"))?,
        artifact_id: match opt_str(map, "artifact_id") {
            Some(id) => Some(ArtifactId::new(id)?),
            None => None,
        },
        digest: match opt_str(map, "digest") {
            Some(d) => Some(Digest::parse(d)?),
            None => None,
        },
        engine: opt_str(map, "engine").map(str::to_owned),
        engine_version: opt_str(map, "engine_version").map(str::to_owned),
        fallback: opt_bool(map, "fallback", false),
    })
}

impl LandscapeResponse {
    /// Serialize to a `landscape_response.schema.json` document.
    pub fn to_json(&self) -> Result<Value> {
        let mut map = Map::new();
        map.insert("schema_version".into(), json!(SCHEMA_VERSION));
        map.insert("contract_major".into(), json!(CONTRACT_MAJOR));
        map.insert("query_id".into(), json!(self.query_id));
        map.insert("status".into(), json!(self.status.as_str()));

        if let Some(problem) = &self.problem {
            map.insert("problem".into(), problem.clone());
        }

        if !self.terms.is_empty() {
            let mut terms = Map::new();
            for (term, value) in &self.terms {
                let mut entry = Map::new();
                entry.insert("value".into(), json!(value.value));
                entry.insert("unit".into(), json!(value.unit.as_str()));
                entry.insert("direction".into(), json!(value.direction.as_str()));
                entry.insert("capability".into(), json!(value.capability.as_str()));
                entry.insert("provenance".into(), provenance_to_json(&value.provenance));
                if let Some(u) = &value.uncertainty {
                    entry.insert("uncertainty".into(), uncertainty_to_json(u));
                }
                if let Some((lo, hi)) = value.reference_range {
                    entry.insert("reference_range".into(), json!([lo, hi]));
                }
                terms.insert(term.as_str().to_owned(), Value::Object(entry));
            }
            map.insert("terms".into(), Value::Object(terms));
        }

        if let Some(vector) = &self.objective_vector {
            map.insert(
                "objective_vector".into(),
                json!({
                    "order": vector.order().iter().map(|t| t.as_str()).collect::<Vec<_>>(),
                    "values": vector.values(),
                    "sense": vector.sense().iter().map(|d| d.as_str()).collect::<Vec<_>>(),
                }),
            );
        }

        if let Some(deltas) = &self.deltas {
            map.insert(
                "deltas".into(),
                json!({
                    "base_digest": deltas.base_digest.as_str(),
                    "exact": deltas.exact,
                    "entries": deltas.entries.iter().map(|entry| json!({
                        "codon_index": entry.codon_index,
                        "alternatives": entry.alternatives.iter().map(|alt| {
                            let mut a = Map::new();
                            if let Some(codon) = alt.codon {
                                a.insert("codon".into(), json!(codon.to_string_upper()));
                            }
                            let mut deltas = Map::new();
                            for (term, d) in &alt.delta {
                                deltas.insert(term.as_str().to_owned(), json!(d));
                            }
                            a.insert("delta".into(), Value::Object(deltas));
                            Value::Object(a)
                        }).collect::<Vec<_>>(),
                    })).collect::<Vec<_>>(),
                }),
            );
        }

        let mut diagnostics = Map::new();
        diagnostics.insert("cache_hits".into(), json!(self.diagnostics.cache_hits));
        diagnostics.insert("cache_misses".into(), json!(self.diagnostics.cache_misses));
        diagnostics.insert("evaluations".into(), json!(self.diagnostics.evaluations));
        if !self.diagnostics.terms_deferred.is_empty() {
            diagnostics.insert(
                "terms_deferred".into(),
                Value::Array(
                    self.diagnostics
                        .terms_deferred
                        .iter()
                        .map(|t| json!(t.as_str()))
                        .collect(),
                ),
            );
        }
        if !self.diagnostics.engine_notes.is_empty() {
            diagnostics.insert("engine_notes".into(), json!(self.diagnostics.engine_notes));
        }
        map.insert("diagnostics".into(), Value::Object(diagnostics));

        Ok(Value::Object(map))
    }

    /// Parse a `landscape_response.schema.json` document.
    ///
    pub fn from_json(value: &Value) -> Result<Self> {
        let map = obj(value, "response")?;
        envelope(map, "response")?;

        if map.contains_key("gradients") {
            return Err(CfError::new(
                code::LAND_CAPABILITY,
                "response carries a gradients block, which no query built by this \
                 contract can request; see ADR-0002",
            ));
        }

        let status = match as_str(map, "status")? {
            "ok" => Status::Ok,
            "partial" => Status::Partial,
            "error" => Status::Error,
            other => return Err(bad(format!("unknown status {other:?}"))),
        };

        let mut terms: Vec<(TermId, TermValue)> = Vec::new();
        if let Some(block) = map.get("terms") {
            for (term, entry) in obj(block, "terms")? {
                let e = obj(entry, "terms[]")?;
                let reference_range = match e.get("reference_range") {
                    None => None,
                    Some(v) => {
                        let a = v
                            .as_array()
                            .ok_or_else(|| bad("reference_range must be an array"))?;
                        if a.len() != 2 {
                            return Err(bad("reference_range must have exactly two entries"));
                        }
                        Some((
                            a[0].as_f64().ok_or_else(|| bad("reference_range[0]"))?,
                            a[1].as_f64().ok_or_else(|| bad("reference_range[1]"))?,
                        ))
                    }
                };
                terms.push((
                    TermId::new(term.as_str())?,
                    TermValue {
                        value: as_f64(e, "value")?,
                        unit: Unit::new(as_str(e, "unit")?),
                        direction: Direction::parse(as_str(e, "direction")?)
                            .ok_or_else(|| bad("unknown direction"))?,
                        capability: match opt_str(e, "capability") {
                            Some(c) => {
                                Capability::parse(c).ok_or_else(|| bad("unknown capability"))?
                            }
                            None => Capability::L0,
                        },
                        provenance: provenance_from_json(get(e, "provenance")?)?,
                        uncertainty: match e.get("uncertainty") {
                            Some(u) => Some(uncertainty_from_json(u)?),
                            None => None,
                        },
                        reference_range,
                    },
                ));
            }
        }

        // A JSON object has no order, so the terms map cannot carry one. The
        // objective vector is what carries the order, and the terms list is
        // sorted here purely so the parse is deterministic.
        terms.sort_by(|a, b| a.0.as_str().cmp(b.0.as_str()));

        let objective_vector = match map.get("objective_vector") {
            None => None,
            Some(v) => {
                let ov = obj(v, "objective_vector")?;
                let order = get(ov, "order")?
                    .as_array()
                    .ok_or_else(|| bad("objective_vector.order must be an array"))?
                    .iter()
                    .map(|t| {
                        let term = TermId::new(
                            t.as_str().ok_or_else(|| bad("objective_vector.order[]"))?,
                        )?;
                        // A term that reached the vector was an axis when the
                        // vector was built; re-asserting that on parse keeps a
                        // hand-edited document from smuggling one in.
                        ObjectiveAxis::new(term, TermRole::ParetoAxis)
                    })
                    .collect::<Result<Vec<_>>>()?;
                let sense = get(ov, "sense")?
                    .as_array()
                    .ok_or_else(|| bad("objective_vector.sense must be an array"))?
                    .iter()
                    .map(|d| {
                        Direction::parse(d.as_str().ok_or_else(|| bad("sense[]"))?)
                            .ok_or_else(|| bad("unknown sense"))
                    })
                    .collect::<Result<Vec<_>>>()?;
                let values = get(ov, "values")?
                    .as_array()
                    .ok_or_else(|| bad("objective_vector.values must be an array"))?
                    .iter()
                    .map(|v| {
                        if v.is_null() {
                            Ok(None)
                        } else {
                            v.as_f64().map(Some).ok_or_else(|| bad("values[]"))
                        }
                    })
                    .collect::<Result<Vec<_>>>()?;
                Some(ObjectiveVector::new(order, sense, values).map_err(bad)?)
            }
        };

        let deltas = match map.get("deltas") {
            None => None,
            Some(v) => {
                let d = obj(v, "deltas")?;
                let entries = get(d, "entries")?
                    .as_array()
                    .ok_or_else(|| bad("deltas.entries must be an array"))?
                    .iter()
                    .map(|entry| {
                        let e = obj(entry, "deltas.entries[]")?;
                        let alternatives = get(e, "alternatives")?
                            .as_array()
                            .ok_or_else(|| bad("alternatives must be an array"))?
                            .iter()
                            .map(|alt| {
                                let a = obj(alt, "alternatives[]")?;
                                let mut delta = Vec::new();
                                for (term, value) in obj(get(a, "delta")?, "delta")? {
                                    delta.push((
                                        TermId::new(term.as_str())?,
                                        value.as_f64().ok_or_else(|| bad("delta value"))?,
                                    ));
                                }
                                delta.sort_by(|x, y| x.0.as_str().cmp(y.0.as_str()));
                                Ok(AlternativeDelta {
                                    codon: match opt_str(a, "codon") {
                                        Some(c) => Some(Codon::parse(c)?),
                                        None => None,
                                    },
                                    delta,
                                })
                            })
                            .collect::<Result<Vec<_>>>()?;
                        Ok(DeltaEntry {
                            codon_index: as_u64(e, "codon_index")? as usize,
                            alternatives,
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;
                Some(Deltas {
                    base_digest: digest(d, "base_digest")?,
                    exact: get(d, "exact")?
                        .as_bool()
                        .ok_or_else(|| bad("deltas.exact must be a boolean"))?,
                    entries,
                })
            }
        };

        let diagnostics = match map.get("diagnostics") {
            None => Diagnostics::default(),
            Some(v) => {
                let d = obj(v, "diagnostics")?;
                Diagnostics {
                    cache_hits: d.get("cache_hits").and_then(Value::as_u64).unwrap_or(0),
                    cache_misses: d.get("cache_misses").and_then(Value::as_u64).unwrap_or(0),
                    evaluations: d.get("evaluations").and_then(Value::as_u64).unwrap_or(0),
                    full_rescores: d.get("full_rescores").and_then(Value::as_u64).unwrap_or(0),
                    rebases: d.get("rebases").and_then(Value::as_u64).unwrap_or(0),
                    terms_deferred: match d.get("terms_deferred") {
                        None => Vec::new(),
                        Some(t) => t
                            .as_array()
                            .ok_or_else(|| bad("terms_deferred must be an array"))?
                            .iter()
                            .map(|t| {
                                TermId::new(t.as_str().ok_or_else(|| bad("terms_deferred[]"))?)
                            })
                            .collect::<Result<Vec<_>>>()?,
                    },
                    engine_notes: match d.get("engine_notes") {
                        None => Vec::new(),
                        Some(n) => n
                            .as_array()
                            .ok_or_else(|| bad("engine_notes must be an array"))?
                            .iter()
                            .map(|n| {
                                n.as_str()
                                    .map(str::to_owned)
                                    .ok_or_else(|| bad("engine_notes[]"))
                            })
                            .collect::<Result<Vec<_>>>()?,
                    },
                }
            }
        };

        Ok(Self {
            query_id: as_str(map, "query_id")?.to_owned(),
            status,
            terms,
            objective_vector,
            deltas,
            diagnostics,
            problem: map.get("problem").cloned(),
        })
    }

    /// RFC 8785 canonical form of the document.
    pub fn to_canonical_json(&self) -> Result<String> {
        canonicalize(&self.to_json()?)
    }
}

impl UncertaintyMethod {
    pub fn parse(text: &str) -> Option<Self> {
        match text {
            "deep_ensemble" => Some(UncertaintyMethod::DeepEnsemble),
            "mc_dropout" => Some(UncertaintyMethod::McDropout),
            "bayesian_last_layer" => Some(UncertaintyMethod::BayesianLastLayer),
            "glm_asymptotic" => Some(UncertaintyMethod::GlmAsymptotic),
            "parametric_bootstrap" => Some(UncertaintyMethod::ParametricBootstrap),
            "split_conformal" => Some(UncertaintyMethod::SplitConformal),
            "mondrian_conformal" => Some(UncertaintyMethod::MondrianConformal),
            "none" => Some(UncertaintyMethod::None),
            _ => None,
        }
    }
}
