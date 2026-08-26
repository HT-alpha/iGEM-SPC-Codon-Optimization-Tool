"""cf-ingest — context-ingestion layer.

FASTA/GenBank parsing to the canonical representation, construct-context assembly
(multi-cistron layout, RBS/promoter annotation), host reference data loaders (tRNA gene
copy number, codon usage, Ribo-seq tracks), and the external-tool seed importers
(LinearDesign / CDSfold / CodonTransformer / DNA Chisel / plain FASTA).

Declared dependencies: cf-types, cf-modelreg.
Contract: planning/02-schemas/construct_context.schema.json
"""
