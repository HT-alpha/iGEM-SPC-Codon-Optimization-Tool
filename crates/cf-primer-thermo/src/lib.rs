//! `cf-primer-thermo` — support layer.
//!
//! Nearest-neighbour Tm and two-state hairpin, self-dimer and hetero-dimer
//! ΔG for short DNA oligos. The parameters are the SantaLucia 1998 Table 2
//! document under `data/`, checked by digest. This crate does not call
//! ViennaRNA and does not depend on the codon search.
//!
//! Declared dependencies: cf-types.
//! Contract: `planning/02-schemas/primer_parameters.schema.json`.

#![forbid(unsafe_code)]

pub mod duplex;
pub mod parameters;
pub mod structure;

pub use duplex::{
    dg_zone, gc_clamp, gc_percent, longest_homopolymer, perfect_duplex, DgZone, Duplex,
};
pub use parameters::{bundled, parse, Parameters};
pub use structure::{dimer, hairpin, Dimer, Fold};

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> Parameters {
        bundled().expect("bundled SantaLucia 1998 document")
    }

    #[test]
    fn table2_cgttga_matches_the_hand_sum() {
        // SantaLucia 1998 Fig. 1 duplex 5'-CGTTGA-3' / 3'-GCAACT-5'.
        // Stacks CG/GC, GT/CA, AA/TT, CA/GT, GA/CT plus one G·C and one A·T initiation.
        // ΔH = -10.6-8.4-7.9-8.5-8.2+0.1+2.3 = -41.2 kcal/mol
        // ΔS = -27.2-22.4-22.2-22.7-22.2-2.8+4.1 = -115.4 cal/K·mol
        let duplex = perfect_duplex(&params(), "CGTTGA").unwrap();
        assert!((duplex.delta_h_kcal - (-41.2)).abs() < 1e-9);
        assert!((duplex.delta_s_cal - (-115.4)).abs() < 1e-9);
        let hand_g = -41.2 - 310.15 * (-115.4 / 1000.0);
        assert!((duplex.delta_g37_kcal - hand_g).abs() < 1e-9);
        // Table 1's rounded ΔG°37 sum is -5.35. Table 2 is coarser (0.1 kcal
        // and 0.1 cal), so the Eq. 2 value differs by about 0.06 kcal/mol.
        assert!((duplex.delta_g37_kcal - (-5.35)).abs() < 0.5);

        let ct = 1e-6_f64 / 4.0;
        let tm_k = (-41.2 * 1000.0) / (-115.4 + 1.987 * ct.ln());
        assert!((duplex.tm_celsius - (tm_k - 273.15)).abs() < 1e-6);
    }

    #[test]
    fn a_replaced_parameter_document_is_refused() {
        let mut text = include_str!("../data/santalucia1998.json").to_string();
        text = text.replace("-7.9", "-7.8");
        let error = parse(&text).expect_err("digest must fail");
        assert!(error.0.contains("CF-PRIM-005"), "{error}");
    }

    #[test]
    fn dg_zones_match_the_idt_cutoffs() {
        assert_eq!(dg_zone(-9.01), DgZone::Refuse);
        assert_eq!(dg_zone(-9.0), DgZone::Flag);
        assert_eq!(dg_zone(-6.0), DgZone::Flag);
        assert_eq!(dg_zone(-4.0), DgZone::Accept);
        assert_eq!(dg_zone(-3.9), DgZone::Accept);
    }

    #[test]
    fn a_gc_run_dimer_is_in_the_refuse_zone() {
        let parameters = params();
        let oligo = "GCGCGCGCGCGCGCGCGCGC";
        let fold = dimer(&parameters, oligo, oligo).unwrap();
        assert!(
            fold.delta_g_kcal < -9.0,
            "GC self-dimer ΔG {} was not below -9",
            fold.delta_g_kcal
        );
        assert_eq!(dg_zone(fold.delta_g_kcal), DgZone::Refuse);
        assert!(fold.alignment.contains('|'));
    }

    #[test]
    fn published_oligoanalyzer_sequences_have_a_table2_dimer() {
        // Nurhafitri and Suwirmen, Bioeksperimen Jurnal Penelitian Biologi,
        // 2026. Sequences from their Table 2. OligoAnalyzer hairpin numbers
        // in that table come from mFold (Owczarzy et al. 2008, NAR 36:W163),
        // which uses loop penalties Table 2 does not contain. This test
        // records the Table 2 two-state dimer, not a claim of IDT parity.
        let parameters = params();
        let examples = [
            ("CAGCCACTGATGCCAAGAGA", "CAGCCACTGATGCCAAGAGA"),
            ("TTGTGGAATGGCATGGCTCT", "TTGTGGAATGGCATGGCTCT"),
            ("CAGCCACTGATGCCAAGAGA", "TTGTGGAATGGCATGGCTCT"),
        ];
        for (left, right) in examples {
            let fold = dimer(&parameters, left, right).unwrap();
            assert!(fold.delta_g_kcal.is_finite(), "{}", fold.delta_g_kcal);
        }
        // Forward self-dimer. Two helices (CAG and CTG) in one register,
        // initiation once. Table 2 gives −3.556 kcal/mol, the published
        // OligoAnalyzer cell −3.55 within 0.01.
        let forward = dimer(&parameters, "CAGCCACTGATGCCAAGAGA", "CAGCCACTGATGCCAAGAGA").unwrap();
        assert!(
            (forward.delta_g_kcal - (-3.556)).abs() < 0.01,
            "{}",
            forward.delta_g_kcal
        );
    }

    #[test]
    fn a_second_helix_crosses_the_flag_line() {
        // 5'-CGCAAAGCG-3' with itself. The best single helix is CGC
        // (or GCG): two stacks and initiation at both ends, ΔG°37 = −2.459,
        // which is above the −4 kcal/mol flag line. The register contains
        // both helices and is charged initiation once, so the pair is flagged.
        let parameters = params();
        let single = perfect_duplex(&parameters, "CGC").unwrap();
        assert!((single.delta_g37_kcal - (-2.459)).abs() < 0.01);
        assert_eq!(dg_zone(single.delta_g37_kcal), DgZone::Accept);
        let both = dimer(&parameters, "CGCAAAGCG", "CGCAAAGCG").unwrap();
        assert!(
            both.delta_g_kcal < -4.0 && both.delta_g_kcal >= -9.0,
            "{}",
            both.delta_g_kcal
        );
        assert_eq!(dg_zone(both.delta_g_kcal), DgZone::Flag);
    }

    #[test]
    fn homopolymer_and_gc_are_literal() {
        assert_eq!(longest_homopolymer("AAACAAAAA"), 5);
        assert!((gc_percent("GCGC") - 100.0).abs() < 1e-9);
        assert!(gc_clamp("ATGC"));
        assert!(!gc_clamp("ATGT"));
    }
}
