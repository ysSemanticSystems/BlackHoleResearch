# AGENTS.md — Contributor & AI-Agent Guide

> Read this **before** writing or modifying any code in this repository. The
> rules below apply to humans and AI agents equally. Cursor will load this
> file and the rules in `.cursor/rules/` automatically.

This repository is **scientific software**. The bar is higher than a typical
application: numbers must be honest, formulas must be cited, and every plot
must be reproducible from on-screen metadata. Read `PHASE2_PLAN.md` for the
strategic direction; this file is the operational guide.

---

## 1. Prime directives

1. **Honesty over completeness.** If you cannot calibrate a number, mark it
   as uncalibrated or upper-limit. **Never invent a value.** Returning
   `None` or raising is preferred over returning a wrong number.
2. **Cite every formula.** A new physics function or constant without a
   primary-source reference in its docstring is rejected.
3. **Use `astropy.units` end-to-end.** Bare floats for physical quantities
   are a bug. Conversion happens at boundaries, not in the middle.
4. **Tests precede expansion.** A new feature lands with a regression test.
   "Test will come later" is not a valid PR description.
5. **Streamlit is presentation only.** No physics or I/O in `app.py`.
   `app.py` may only call into `blackhole.*` modules and render.
6. **Single source of truth.** Source parameters live in
   `blackhole/catalog.py`; survey calibrations in `blackhole/calibration.py`;
   physics constants come from `astropy.constants`. **Do not** hardcode
   alternatives elsewhere.
7. **Provenance ships with every plot.** A figure without metadata
   identifying the input file SHA, calibration version, and function chain
   is incomplete.
8. **Pitfalls are normative.** Items in `docs/PITFALLS.md` describe how
   things break in practice. Code that contradicts an entry there is wrong;
   amend the entry only with strong justification.

---

## 2. Repository map

```
BlackHoleResearch/
├── app.py                       # Streamlit UI. Presentation only.
├── pyproject.toml               # Build, deps, tool config (M0)
├── requirements.txt             # Pinned bounds
├── requirements-lock.txt        # Byte-reproducible install (M0)
├── README.md                    # User-facing landing page
├── AGENTS.md                    # This file
├── PHASE2_PLAN.md               # The plan
├── CITATION.cff                 # Machine-readable citation (M11)
├── CHANGELOG.md                 # Keep-a-Changelog format (M0)
├── .cursor/
│   └── rules/                   # Persistent AI guidance, one concern per file
├── .github/
│   └── workflows/ci.yml         # CI matrix (M0)
├── blackhole/
│   ├── __init__.py
│   ├── catalog.py               # Source dataclass + registry           (M1)
│   ├── calibration.py           # Per-survey flux conversion             (M2)
│   ├── photometry.py            # Aperture photometry → SEDPoint         (M3)
│   ├── cosmo.py                 # Planck18-default cosmology helpers     (M7)
│   ├── provenance.py            # Provenance dataclass + helpers         (M6)
│   ├── jwst.py                  # JWST NIRSpec/NIRCam ingestion          (M8)
│   ├── lrd.py                   # Little Red Dot diagnostics              (M9)
│   ├── io.py                    # FITS load/inspect/event-binning
│   ├── wcs_plot.py              # Image rendering with WCS
│   ├── spectra.py               # PHA spectra
│   ├── sed.py                   # SED assembly
│   ├── lightcurves.py           # Time-series binning + LS periodogram
│   ├── _style.py                # Shared matplotlib styling              (M6)
│   ├── _rgb.py                  # Lupton+2004 RGB composites              (M10)
│   └── physics/
│       ├── accretion.py
│       ├── spectral_xray.py
│       ├── infrared.py
│       └── variability.py
├── scripts/
│   ├── download_data.py         # SkyView Phase-1 cutouts
│   └── download_jwst.py         # MAST Phase-2 LRD downloads             (M8)
├── tests/
│   ├── conftest.py              # `tiny_fits`, synthetic data fixtures
│   ├── test_*.py                # One test file per public module
│   └── test_physics/*.py        # One per physics submodule
├── docs/
│   ├── PITFALLS.md              # Failure-mode catalog (normative)
│   ├── adr/                     # Architecture Decision Records
│   ├── calibration_table.md     # Per-survey calibration formulas        (M2)
│   ├── spectral_fitting.md      # Sherpa setup, if M5b is taken
│   ├── literature_seds.py       # Old hardcoded SEDS, kept for overlay   (M3)
│   └── notebooks/
│       ├── 01_ngc1068_seyfert2.ipynb                                     (M11)
│       ├── 02_m87_jet.ipynb                                              (M11)
│       ├── 03_cygx1_xray.ipynb                                           (M11)
│       └── 04_lrd_diagnostics.ipynb                                      (M11)
└── fits_data/
    ├── MANIFEST.json
    └── *.fits
```

Anything not in this map should not exist. If you need a new top-level
location, propose it via an ADR in `docs/adr/`.

---

## 3. How to work on this codebase

### 3.1 Before you start

1. Pull `main`; create a feature branch named after the milestone:
   `m3-aperture-photometry`, `m7-lrd-catalog`, etc.
2. Read the relevant milestone in `PHASE2_PLAN.md` end-to-end. Exit criteria
   are binary.
3. Read the relevant `.cursor/rules/*.mdc` for the file types you'll touch.

### 3.2 Workflow

1. Write the tests first (or alongside) the implementation.
2. Implement against the contract documented in `PHASE2_PLAN.md` §2 and §4.
3. Run locally before pushing:
   ```bash
   ruff check . && ruff format --check .
   mypy --strict blackhole/
   pytest -q --cov=blackhole --cov-fail-under=70
   ```
4. Update `CHANGELOG.md` under `## [Unreleased]`.
5. Open a PR. The PR description must reference the milestone ID and tick
   off the exit criteria.

### 3.3 What gets a PR rejected

- Test coverage drops below the threshold.
- `mypy --strict` regresses.
- A new physics formula without a primary-source citation in the docstring.
- A hardcoded flux, distance, M_BH, or redshift outside `blackhole/catalog.py`
  or `tests/`.
- A new top-level file or directory not in the repo map (§2) without an ADR.
- A `from blackhole.*` removed without a deprecation pass.
- A new dependency added without justification in the PR description.
- Computation in `app.py` outside of trivial value lookups.

---

## 4. Code conventions (summary; see `.cursor/rules/` for detail)

### 4.1 Python

- Python 3.11+ only. Use modern syntax: `X | None`, `list[int]`, `match`.
- `from __future__ import annotations` at the top of every module.
- Dataclasses for data; functions for behavior; classes only when there
  is real state.
- Public dataclasses are `frozen=True` unless mutation is essential.
- Public functions are fully type-annotated, including return types.
- Private helpers begin with `_`.
- Module-level mutable state is forbidden except for caches and registries.

### 4.2 Imports

```python
# Standard library
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path

# Third-party scientific
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS

# First-party
from blackhole.catalog import Source, by_short_id
from blackhole.io import ImageData
```

Three groups, alphabetised within each, single blank lines between.

### 4.3 Units

Anything that has a physical dimension carries a `Quantity`:

```python
# GOOD
flux_density: u.Quantity = 1.6 * u.Jy
wavelength: u.Quantity = 2.2 * u.micron

# BAD — drop on the floor
flux_density: float = 1.6
```

Conversion happens at:
- The boundary with the FITS file (calibrator-side).
- The plotting axis (renderer-side).
- Inputs to scipy/numpy that don't accept Quantity (use `.to(u.X).value`).

Never inside physics.

### 4.4 Errors

- Calibration on uncalibratable data → `UncalibratedDataError`.
- Missing required header keyword → `FITSHeaderError`, with the missing
  keyword named.
- Empty selection (e.g., zero events after filter) → return an empty
  container, **not** raise.

### 4.5 Logging

- Use Python `logging` at the module level, never `print()` (except
  in `scripts/`).
- Logger name = module name: `logger = logging.getLogger(__name__)`.

### 4.6 Docstrings

The shape every module-level docstring must have:

```python
"""
blackhole.<module> — one-line summary.

WHY THIS MODULE EXISTS
======================
2–6 lines of context. What problem does this solve? Why is it a module
instead of a few lines elsewhere?

KEY PHYSICS / CONVENTIONS
=========================
3–10 lines. Define the canonical quantities and equations this module
implements. Cite where the convention comes from.

PITFALLS
========
Bullet list of items from docs/PITFALLS.md that apply here.

REFERENCES
==========
- Primary source 1 (Author Year, journal, vol, page or DOI/arXiv).
- Primary source 2 ...
"""
```

The shape every public-function docstring must have:

```python
def fn(arg: T) -> R:
    """One-line summary.

    Longer prose where useful.

    Parameters
    ----------
    arg : T
        What it is.

    Returns
    -------
    R
        What you get.

    References
    ----------
    Author Year, journal, vol, page (or DOI/arXiv).
    """
```

---

## 5. The honesty checklist (apply to every PR that produces numbers)

Before merging a PR that adds or changes a measurement function or a UI
display of a number:

- [ ] **Units.** Every quantity displayed has explicit units.
- [ ] **Uncertainty.** Every measurement has a 1-σ error or a stated reason
      why one is not available (and that reason is shown to the user).
- [ ] **Upper limits.** If the measurement is non-detection (S/N < 3),
      it's marked as such.
- [ ] **Provenance.** The plot or table includes the source file SHA, the
      calibration version, and the function chain.
- [ ] **Reference.** Any threshold or formula cites a primary source.
- [ ] **Test.** A unit test exercises the calculation with a hand-computed
      expected value.

---

## 6. Git hygiene

- Commit messages: imperative mood, subject < 72 chars, body wrapped at 80.
- One logical change per commit; the PR may contain multiple commits.
- Squash on merge unless commits tell a useful story.
- Never commit `.DS_Store`, `__pycache__/`, `.joblib_cache/`, `outputs/`,
  or `fits_data/*.fits` (already in `.gitignore`).
- Never commit credentials, API keys, or `.env` files.

---

## 7. When in doubt

1. Read `PHASE2_PLAN.md` — the answer to most "should I do X?" is in §4.
2. Read `docs/PITFALLS.md` — the answer to most "is X a known problem?" is
   in §1–14.
3. Read the relevant rule in `.cursor/rules/` — the answer to most
   "how should I format X?" is there.
4. Write an ADR in `docs/adr/` proposing the change.

If none of those help, leave a `# TODO(name): reason and ticket` and surface
it in the PR description.
