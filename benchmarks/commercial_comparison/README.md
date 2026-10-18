# DNA Chisel comparison runner

Reproduces the runnable arm of
[`docs/benchmarks/commercial-codon-tools.md`](../../docs/benchmarks/commercial-codon-tools.md).

Commercial web tools are not called. The gene list and the job settings are
fixed in `run_chisel_vs_dds.py`. They match the study in
`docs/benchmarks/v2.9.3-external-comparison.md`.

From the repository root, after a release build:

```sh
.venv/bin/maturin develop --release
.venv/bin/pip install 'dnachisel==3.2.16' 'ViennaRNA==2.7.2'
.venv/bin/python benchmarks/commercial_comparison/run_chisel_vs_dds.py \
    --out docs/benchmarks/commercial-comparison-rows.json
```

The script writes the JSON rows and
`docs/benchmarks/commercial-comparison-tables.md`.
