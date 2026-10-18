"""cf-biophys-lm — biophysical-predictors layer.

Codon language-model naturalness term (CodonTransformer / CodonBERT class models).
L0-only by default, optional L2 gradients w.r.t. relaxed codon logits. Runs in the batch
worker because it needs torch and optionally a GPU.

Declared dependencies: cf-types, cf-landscape-api, cf-modelreg.
Contract: planning/02-schemas/interfaces/predictors.pyi
"""
