"""DNA Design Suite — codon optimization and primer design.

This package is the Python half of the two-runtime split decided in
`planning/07-decisions/ADR-0001-language-and-runtime-split.md`: model training,
data ingestion, job orchestration and the API surface live here, while the
per-move hot path lives in the compiled core exposed as ``dna_design_suite._core``.

The boundary rule from ADR-0001 — no call crosses it more than once per
generation — is enforced mechanically by ``tools/check_module_layering.py``,
which confines the pyo3 binding surface to a single crate.
"""

from __future__ import annotations

__version__ = "0.1.0.post1"

__all__ = ["__version__", "core_version", "core_available"]


def core_available() -> bool:
    """Whether the compiled core extension is importable.

    Pure-Python parts of the package (schema handling, manifests, ingest) are
    usable without it, which keeps ``pip install -e .`` working for a
    contributor who has not built the Rust core yet.
    """
    try:
        import dna_design_suite._core  # noqa: F401
    except ImportError:
        return False
    return True


def core_version() -> str:
    """Version string reported by the compiled core.

    Raises ``ImportError`` with an actionable message rather than a bare
    ``ModuleNotFoundError`` when the core has not been built.
    """
    try:
        from dna_design_suite import _core
    except ImportError as exc:  # pragma: no cover - exercised in CI without core
        raise ImportError(
            "dna_design_suite._core is not built. Run `maturin develop -m "
            "crates/cf-pyext/Cargo.toml` in a development checkout, or install "
            "a released wheel, which ships the core prebuilt."
        ) from exc
    return _core.version()


def design(*args, **kwargs):
    """`cf.design(...)` — the documented library entry point.

    Imported lazily so that `import dna_design_suite` stays cheap and so that the
    top-level package does not pull the whole orchestration stack in for a
    caller who only wanted `_core`.
    """
    from dna_design_suite.codon.orchestrator import design as _design

    return _design(*args, **kwargs)
