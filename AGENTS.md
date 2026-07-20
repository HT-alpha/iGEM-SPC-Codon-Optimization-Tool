# AGENTS.md

## Cursor Cloud specific instructions

### Current repository state

As of this writing, this repository is a **stub**: the only tracked file is `README.md`
(a single-line title, `# iGEM Codon Optimization Tool`). There is:

- No application/source code
- No dependency manifest (`package.json`, `requirements.txt`, `pyproject.toml`, lockfiles, etc.)
- No build, test, lint, or run tooling
- No services, databases, or external dependencies

There is therefore nothing to install, build, run, or test yet. The startup update
script is intentionally a no-op until real project tooling is added.

### Available runtimes in the Cloud VM

The base image already provides (no install needed): Node.js 22.x + npm 10.x, and
Python 3.12 + pip. Use these when bootstrapping the actual project.

### When application code is added

Once real project files land, revisit environment setup:

- Add the dependency-install step to the startup update script (e.g. `npm ci` /
  `pip install -r requirements.txt` / `uv sync`), matching whatever package manager
  the lockfile implies.
- Document the dev/run/lint/test commands here (or reference `package.json` scripts /
  `Makefile` targets rather than duplicating them).
