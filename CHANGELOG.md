# Changelog

All notable changes to this project are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added — M2 calibration registry
- `blackhole/calibration.py` — frozen `CalibratedImage` dataclass,
  `UncalibratedDataError`, `CALIBRATION_VERSION = "1.0.0"`,
  header-first `detect_survey` dispatch, and `calibrate(image)`.
- Per-survey calibrators: 2MASS-K (Cohen+2003 F_nu_0=666.7 Jy),
  IRIS_12/25/60/100 (MJy/sr → Jy/pixel via solid angle), AKARI
  (Doi+2015), VLA-FIRST (Jy/beam → Jy/pixel via BMAJ/BMIN beam area,
  Becker+1995), RASS (counts/EXPTIME * ECF=1.08e-11 from Snowden+1995),
  and DSS which raises `UncalibratedDataError`.
- `tests/test_calibration.py` — 18 tests, per-survey 2-pixel synthetic
  FITS matched to hand-computed Quantity within 1e-4.
- `docs/calibration_table.md` — per-survey conversion table, references,
  trap notes, and calibration-version bump policy.

### Added — M3 aperture photometry pipeline
- `blackhole/photometry.py` — `PhotometryResult` frozen dataclass and
  `aperture_photometry_on(cal, coord, ...)`. WCS-aware sky-aperture
  placement, sigma-clipped MAD sky stats on the annulus, 3-sigma
  upper-limit rule.
- `aperture_for_band(band)` — sensible defaults (radio 60", submm 120",
  IR 20", opt/UV 10", X-ray 30").
- `tests/test_photometry.py` — 14 tests: 2D Gaussian recovery within 2%
  through the 2MASS-K calibration, upper-limit detection, Jy and
  erg/s/cm² routing in `to_sed_point`, WCS-less fallback.
- `docs/literature_seds.py` — literature SED and X-ray reference values
  moved out of `app.py` and packaged with full bibcodes.
- **app.py** SED tab — local cutouts photometred via the new pipeline
  (opt-in via "Aperture-photometer the local FITS cutouts"); literature
  values are an opt-in overlay tagged `[lit]`. Local-photometry
  diagnostics expander surfaces per-file flux, error, aperture, and
  skip reasons (DSS, etc).
- `pyproject.toml` — `photutils>=2.0` added to runtime dependencies.

### Added — M4 GTI-aware light curves
- `blackhole/io.py` — `EventList.gti` field (N×2 start/stop seconds);
  `_read_gti(hdul)` reads GTI, STDGTI, and per-CCD/per-FPM `GTI*`
  extensions and unions overlapping intervals.
- `blackhole/lightcurves.py` — `LightCurve` carries `effective_exposure`,
  `gti_applied`, `total_exposure_s`, `gti_total_s`.
  `bin_events_to_lightcurve(..., apply_gti=True, min_exposure_fraction=0.5)`
  spans the full GTI window, drops bins below the threshold, and
  normalizes rate and Poisson error by effective exposure
  (Vaughan+2003).
- `render_lightcurve` title now reports GTI exposure and fraction.
- `tests/conftest.py` — `gapped_events_fits` fixture (two 1 ks GTIs with
  a 1 ks gap).
- `tests/test_lightcurves.py` — 9 GTI tests: no fake zero bins in the
  gap, ≈ 2000 s effective exposure, partial-bin rate normalization,
  title metadata, and a non-GTI fallback path.

### Changed — M5a spectrum tab honesty
- `blackhole/spectra.py` — `fit_power_law` returns α_channel (descriptive
  channel-space slope), not Γ. `render_spectrum` x-axis is
  "Channel (not energy)"; legend label is `α_channel`.
- **app.py** Spectrum tab — `classify_photon_index(g)` call *removed*
  from the channel-space code path. Subheader retitled
  "1D X-ray spectrum (channel space)". An expander explains why this
  is α_channel and not Γ, with a pointer to M5b.
- `tests/test_spectra.py` — 3 new tests pin the contract: the xlabel
  must contain both "Channel" and "not energy"; the legend says
  `α_channel`; `classify_photon_index` is patched to raise, and the
  channel-space path still succeeds — proving it never touches the
  classifier.

### Added — M6 provenance & UI revamp
- `blackhole/provenance.py` — frozen `Provenance` dataclass (fits_sha256,
  fits_path, calibration_version, function_chain, library_version,
  timestamp_utc, schema_version, extra), JSON round-trip,
  `sha256_of_file` (returns "synthetic" for missing files),
  `build_provenance`, `attach`/`get`/`extend_chain`, `save_figure`
  (PNG + sidecar JSON), `load_sidecar`, `as_table_rows`,
  `as_bibtex_note`.
- `blackhole/_style.py` — `DARK_BG`, `DARK_FG`, `apply_dark(ax)`,
  `apply_dark_figure(fig)`. Single source of truth for the dark-theme
  spine / tick / grid / label styling.
- Every `render_*` in `blackhole.*` now routes styling through
  `apply_dark` and attaches a `Provenance` to its returned Figure.
- **app.py** — `show_with_provenance(fig, key, copy_target)` helper used
  by every tab: renders the figure, surfaces a Provenance table plus a
  copy-citation snippet, and optionally writes a PNG + sidecar JSON to
  `outputs/` when the sidebar toggle is on.
- `tests/test_provenance.py` — 15 tests: schema round-trip, SHA-256
  for real and missing files, every renderer attaches a Provenance
  with the right function_chain, save_figure round-trip, end-to-end
  calibrated photometry → SED render carries `calibration_version`,
  and a grep test enforces that the inline
  `for spine in ax.spines.values(): spine.set_color(...)` idiom lives
  only inside `_style.apply_dark`.

### Added — M1 source catalog
- `blackhole/catalog.py` — frozen `Source` dataclass + immutable
  `CATALOG` tuple for NGC 1068, M87, Cyg X-1. Every numeric field
  carries a primary-source reference string.
- Helpers: `by_short_id`, `by_name`, `by_filename`,
  `eddington_luminosity_of`, `distance_to`, `redshift_of`.
- `tests/test_catalog.py` — 42 tests pinning each value to its
  citation: NGC 1068 M_BH to Lodato & Bertin 2003, M87 to EHT 2019
  Paper VI, Cyg X-1 to Miller-Jones+2021. Eddington luminosities
  hand-checked against Frank+2002 Eq. 1.5.
- **app.py** — catalog-driven science banner under the file metric
  strip (Source / Type / M_BH / Distance / L_Edd) plus a per-source
  References expander. SED tab's target list and the filename-to-target
  resolver now flow from `blackhole.catalog` (removed the local
  `TARGET_FILENAME_KEYS` dict and the hardcoded `["NGC 1068", ...]`
  array).

### Changed — UX: Overview tab and per-tab display controls
- Added an **Overview** tab as the new landing surface: a 3-column
  thumbnail grid that renders every file in `fits_data/` (image HDUs
  via stretched intensity, event-list files via low-resolution bin maps).
  Cached per (path, stretch, cmap).
- Moved the **file picker** out of the sidebar and into a top-of-page
  selector next to the metric strip, so file choice and file metadata
  sit together.
- Moved **stretch and colormap** controls out of the sidebar and into
  per-tab "Display options" expanders on the Image and X-ray Events
  tabs. Each tab now owns the display state it actually drives.
- The SED tab now **auto-detects** its target from the active filename
  (recognises `ngc1068`, `m87`, `cygx1` patterns), preselecting the
  matching SED entry so the page lands ready to read.
- Sidebar trimmed to **orientation content only**: dataset stats,
  target list, data-source links, and repo/plan/pitfalls cross-links.

### Deferred — explicitly tracked for M1
- "Science banner" redesign (target-aware physics summary above the
  fold). Currently kept as a generic intro block; will become
  catalog-driven once the source-catalog work in M1 lands. Tracked in
  `PHASE2_PLAN.md` (M1 scope notes) and as a comment at the top of
  `app.py`.

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
