# Developer entry points. CI runs the same commands (.github/workflows/ci.yml),
# so a green `make check` locally means a green CI run.

VENV ?= .venv
PY   := $(VENV)/bin/python
PIP  := $(VENV)/bin/pip

.PHONY: release-check help setup build test check fmt lint invariants rust python assets clean

help:
	@echo "setup      create the venv and install dev dependencies"
	@echo "build      compile the Rust core into the venv (maturin develop)"
	@echo "test       cargo test + pytest"
	@echo "check      everything CI runs"
	@echo "invariants planning package, module layering, schema sync"
	@echo "assets     regenerate vendored assets and show any diff"

# `python3 -m venv` needs `ensurepip`, which Debian and Ubuntu ship in a
# separate `python3-venv` package that is not always present and is not always
# installable. `virtualenv` vendors its own bootstrap and needs neither, so it
# is tried first and installed if missing; the stdlib module remains the last
# resort rather than the first. Getting this wrong looks like a broken clone
# rather than a missing package, which is the wrong first impression for a
# contributor following the README.
setup:
	python3 -m virtualenv -q $(VENV) \
	  || (python3 -m pip install --quiet --user --break-system-packages virtualenv \
	      && python3 -m virtualenv -q $(VENV)) \
	  || python3 -m venv $(VENV)
	$(PIP) install -q --upgrade pip
	$(PIP) install -q "maturin>=1.5,<2" pytest jsonschema pyyaml biopython ruff

build:
	$(VENV)/bin/maturin develop

rust:
	cargo test --workspace --all-targets

python: build
	$(PY) -m pytest -q

test: rust python

fmt:
	cargo fmt --all
	$(VENV)/bin/ruff format python tools planning/tools

lint:
	cargo fmt --all -- --check
	cargo clippy --workspace --all-targets -- -D warnings
	$(VENV)/bin/ruff check python tools planning/tools

invariants:
	$(PY) planning/tools/validate_planning_package.py
	$(PY) tools/check_module_layering.py --check
	$(PY) tools/check_schema_sync.py --check

assets:
	$(PY) tools/generate_genetic_codes.py
	@git diff --stat -- crates/cf-types/data/genetic_codes.json

release-check:
	@# AMD-013: the same gate CI runs before building a release wheel.
	$(PY) tools/check_release_clearance.py

check: invariants lint test

clean:
	cargo clean
	rm -rf $(VENV) .pytest_cache
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
