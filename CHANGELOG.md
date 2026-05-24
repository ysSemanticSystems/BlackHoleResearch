# Changelog

All notable changes to this project are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added — M0 regression net & repo hardening
- `pyproject.toml` — PEP 621 project metadata, dependency pins, pytest and
  coverage config; installable via `pip install -e ".[dev]"`.
- `ruff.toml` — lint configuration (`E,W,F,I,B,UP,SIM,RUF,N`) with explicit
  carve-outs for physics symbol naming and Greek/Unicode scientific text.
- `mypy.ini` — strict type-checking config; library-side third-party stubs
  declared per-module.
- `tests/conftest.py` — `tiny_image_fits`, `tiny_events_fits`, `tiny_pha_fits`
  fixtures: synthetic, deterministic FITS files for offline test runs.
- Seed tests: `tests/test_io.py`, `tests/test_sed.py`, `tests/test_spectra.py`,
  `tests/test_lightcurves.py`, `tests/test_wcs_plot.py`, and four physics
  test modules under `tests/test_physics/`. **100 tests total, 92% line
  coverage** (gate at 70%).
- `.github/workflows/ci.yml` — matrix CI: pytest on Linux+macOS × Python
  3.11+3.12, plus dedicated ruff and mypy jobs; coverage uploaded as an
  artifact on the Linux/py3.12 leg.

### Fixed
- `blackhole/io.load_image`: previously accepted a header with no celestial
  CTYPE keywords as a valid WCS (astropy silently returns an identity 2-axis
  WCS). Now requires `WCS.has_celestial` before returning a non-None WCS;
  matches the rule in `docs/PITFALLS.md` #3.

### Added — Discoverability & metadata
- `LICENSE` — explicit MIT license file (README previously claimed MIT
  without one).
- `CITATION.cff` — Citation File Format v1.2.0 with full keyword set,
  references to primary-source papers, and GitHub citation-widget support.
- `codemeta.json` — CodeMeta 2.0 research-software metadata (JSON-LD).
- `.zenodo.json` — Zenodo metadata for the planned v0.2.0 DOI mint.
- `README.md` — comprehensive overhaul for researcher discoverability:
  badges row, value-prop above the fold, "why you'll want this" entry
  points, search-keyword block (mission, technique, formula names),
  feature deep-dives with primary-source citations, related-tool index
  (DS9, CIAO, XSPEC, SAS, Stingray, aplpy), standards-adhered-to list,
  BibTeX citation block, acknowledgements.

### Added — Phase 2 governance
- `PHASE2_PLAN.md` — milestone roadmap with binary exit criteria, no time
  estimates, dependency graph, and per-milestone scope.
- `AGENTS.md` — operational guide for contributors and AI agents.
- `.cursor/rules/` — single-concern persistent rules:
  - `00-project-context.mdc` (always loaded)
  - `01-documentation-standards.mdc`
  - `02-scientific-honesty.mdc`
  - `03-python-style.mdc`
  - `04-fits-handling.mdc`
  - `05-testing-discipline.mdc`
  - `06-ui-patterns.mdc`
  - `07-reproducibility.mdc`
  - `08-adr.mdc`
  - `09-commit-and-pr.mdc` (always loaded)
- `docs/adr/` — Architecture Decision Records directory with seed index.
- `CHANGELOG.md` (this file).

### Notes
- The Phase 2 milestone work has not yet started; the items above are
  governance and planning scaffolding only.

---

## [0.1.0] — Phase 1 (2026-05)

### Added
- Streamlit multiwavelength FITS explorer (`app.py`).
- `blackhole/` library: `io`, `wcs_plot`, `spectra`, `sed`, `lightcurves`.
- `blackhole/physics/`: `accretion`, `spectral_xray`, `infrared`,
  `variability`.
- `scripts/download_data.py` — SkyView cutouts for NGC 1068, M87, Cyg X-1.
- `docs/PITFALLS.md` — failure-mode catalog (14 entries).
- `fits_data/MANIFEST.json` — SHA-256-stamped manifest of downloads.
