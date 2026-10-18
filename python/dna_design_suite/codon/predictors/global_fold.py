"""Optional whole-sequence fold, labelled apart from windowed accessibility.

`AMD-001` keeps global minimum free energy out of the per-move search. This
module is not on `StructureEngine`. Nothing in the tier loop calls it.

When ViennaRNA is installed, the value is `RNA.fold`'s minimum free energy
in kcal/mol, Turner nearest-neighbour parameters as that install ships them,
scope `global_mfe`. It is not `structure.local_accessibility`.

When ViennaRNA is not installed, a sequence of at most 400 nt gets a
Nussinov maximum matching: Watson-Crick pairs only, hairpin loop of at
least 3 unpaired bases, no stacking energy. The value is a base-pair count.
It is not a free energy and not a Turner ΔG. Longer sequences are not
scored with that cubic recursion, and the report says so instead of
inventing a number.

`DEP-STRU-007`.
"""

from __future__ import annotations

from typing import Any

_NUSSINOV_MAX_NT = 400
_PAIRS = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}


def global_fold_report(sequence: str, *, temperature_c: float = 37.0) -> dict[str, Any]:
    """One labelled global fold of `sequence`."""
    dna = "".join(sequence.upper().split())
    if not dna or any(base not in "ACGTU" for base in dna):
        raise ValueError("global fold needs an A/C/G/T sequence")
    rna = dna.replace("T", "U")
    report: dict[str, Any] = {
        "schema_version": "1.0.0",
        "kind": "global_fold",
        "length_nt": len(rna),
        "scope": "whole_sequence",
        "not_windowed_accessibility": True,
        "not_a_vendor_synthesis_score": True,
        "wet_lab_expression": "out_of_scope",
        "temperature_c": temperature_c,
    }
    vienna = _vienna_mfe(rna, temperature_c)
    if vienna is not None:
        structure, mfe = vienna
        report.update(
            {
                "engine": "viennarna",
                "quantity": "global_mfe_kcal_per_mol",
                "unit": "kcal/mol",
                "value": mfe,
                "dot_bracket": structure,
                "label": (
                    "ViennaRNA global minimum free energy (kcal/mol) of the whole sequence. "
                    "This is not windowed unpaired probability, and it is not "
                    "structure.local_accessibility."
                ),
                "fallback": None,
            }
        )
        return report
    fallback = _nussinov_or_refusal(rna)
    report.update(
        {
            "engine": fallback["engine"],
            "quantity": fallback["quantity"],
            "unit": fallback["unit"],
            "value": fallback["value"],
            "dot_bracket": None,
            "label": fallback["label"],
            "fallback": fallback,
            "viennarna": "not_installed",
        }
    )
    return report


def _vienna_mfe(rna: str, temperature_c: float) -> tuple[str, float] | None:
    try:
        import RNA
    except ImportError:
        return None
    from dna_design_suite.codon.predictors.viennarna import FOLD_LOCK

    with FOLD_LOCK:
        RNA.cvar.temperature = float(temperature_c)
        structure, mfe = RNA.fold(rna)
    return structure, float(mfe)


def _nussinov_or_refusal(rna: str) -> dict[str, Any]:
    if len(rna) > _NUSSINOV_MAX_NT:
        return {
            "engine": "unavailable",
            "quantity": None,
            "unit": None,
            "value": None,
            "label": (
                f"ViennaRNA is not installed, and the Nussinov fallback is only "
                f"computed up to {_NUSSINOV_MAX_NT} nt because it is cubic. "
                "No free energy was invented."
            ),
        }
    pairs = nussinov_pair_count(rna)
    return {
        "engine": "nussinov_fallback",
        "quantity": "nussinov_watson_crick_pairs",
        "unit": "base_pairs",
        "value": pairs,
        "label": (
            "Nussinov maximum Watson-Crick matching (A-U and G-C only, hairpin "
            "loop at least 3). This is a base-pair count, not a kcal/mol minimum "
            "free energy, and not a Turner parameter set. ViennaRNA is not installed."
        ),
    }


def nussinov_pair_count(rna: str) -> int:
    """Maximum number of Watson-Crick pairs. No energy parameters."""
    seq = rna.upper().replace("T", "U")
    n = len(seq)
    if n == 0:
        return 0
    dp = [[0] * n for _ in range(n)]
    for span in range(1, n):
        for i in range(n - span):
            j = i + span
            best = dp[i][j - 1]
            if j - i > 3 and (seq[i], seq[j]) in _PAIRS:
                inner = dp[i + 1][j - 1] if i + 1 <= j - 1 else 0
                best = max(best, inner + 1)
            for k in range(i, j):
                right = dp[k + 1][j] if k + 1 <= j else 0
                best = max(best, dp[i][k] + right)
            dp[i][j] = best
    return dp[0][n - 1]
