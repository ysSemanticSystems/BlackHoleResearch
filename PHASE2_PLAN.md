# BlackHoleResearch — Phase 2 Plan

> A milestone-based roadmap to turn this repository from a *teaching demo* into a
> *citable multiwavelength research instrument*, then extend it to the Little Red
> Dot frontier. **No time estimates.** Each milestone is gated by a testable
> exit criterion and an explicit dependency graph. Work is finished when the
> exit criteria pass on `main`, not when "it feels done."

---

## 0. North Star

### 0.1 Vision
A reproducible, well-cited, calibration-aware multiwavelength explorer that:
1. Loads real FITS data from heterogeneous missions.
2. Produces **calibrated** measurements with **stated uncertainties** and
   **traceable provenance** for every number rendered to the user.
3. Hosts a focused Phase-2 module for **Little Red Dots (LRDs)** that lets a
   researcher reproduce the canonical "LRD violates AGN unification" plot
   (L_X/L_UV vs Lusso & Risaliti 2016) from raw archival files.
4. Ships with a DOI, a regression test suite, and a worked end-to-end notebook
   per target.

### 0.2 What "great" means — operationally
A neutral observer can answer **yes** to all of:

- [ ] Can I `git clone`, `pip install -e .`, `python scripts/download_data.py`,
      `pytest`, and `streamlit run app.py` on a clean machine without manual
      patches?
- [ ] If I delete the hardcoded SED lookup table, do the SED tab numbers come
      from aperture photometry on the local FITS cutouts?
- [ ] Is every flux value rendered in the UI accompanied by (a) units,
      (b) a 1-σ uncertainty, (c) a provenance string (file SHA + calibration
      version + function name)?
- [ ] Does the test suite cover ≥80% of `blackhole/` and ≥60% of the
      photometry/calibration paths, with at least one *physics regression*
      test (e.g. Eddington luminosity for 1 M☉ → 1.26e38 erg s⁻¹ to 4 sig
      figs)?
- [ ] Does CI run on Python 3.11 and 3.12 with `ruff`, `mypy --strict`, and
      `pytest` green on every PR?
- [ ] Is there a Zenodo DOI on a tagged release?
- [ ] Can I reproduce **one published L_X/L_UV value** for an LRD from
      Greene+2024, Maiolino+2024, or Akins+2025 within the stated literature
      tolerance, using the notebook in `docs/notebooks/`?

### 0.3 Non-goals (explicit)
- We will **not** become a competitor to XSPEC, Sherpa, or CIAO. Their fits are
  the gold standard; we wrap or link them where it matters and stay out of
  their lane otherwise.
- We will **not** ingest raw level-1 mission data (event reprocessing, flagged
  pixel masks, dark current). We consume **standard archival products**.
- We will **not** attempt to image-model jets at the EHT/VLBI level. That is
  Phase 3 (or someone else's project).
- We will **not** support Python < 3.11.

### 0.4 Audience
Three concentric rings:
1. **Self** — internal research instrument.
2. **Students and observatories** — teaching tool with reproducible labs.
3. **Reviewers and collaborators** — referenced in a methods section as the
   tool used to produce a figure.

The closer to ring 3 we push, the more the "hardening" pillar matters.

### 0.5 Citable-artifact criterion
A reference of the form:

> SEDs and L_X/L_UV ratios were computed using **BlackHoleResearch v0.2.0**
> (Ye, 2026), DOI 10.5281/zenodo.XXXXXXX.

must be defensible. That is the bar.

---

## 1. Strategic principles

These are the load-bearing decisions; every milestone is consistent with them.

| # | Principle | Consequence in code |
|---|---|---|
| P1 | **Honesty over completeness.** A missing or upper-limit value is a feature, not a bug. | `SEDPoint.upper_limit` exists; we propagate `None` through `astropy.units.Quantity` chains; we never invent a number. |
| P2 | **Close the data ↔ analysis loop.** Hardcoded numbers in the UI are a tech debt to be paid down, not a permanent feature. | A test asserts that no module under `blackhole/` contains a literal flux in Jy outside of `tests/` and `docs/`. |
| P3 | **Single source of truth for sources.** Every target lives in `blackhole.catalog`, with M_BH, distance, redshift, type, references. | UI and physics both import from the catalog; the catalog is unit-tested against published values. |
| P4 | **Calibration is a first-class concern.** No "arbitrary units" leaves a measurement function. | Every photometry function returns either a calibrated `Quantity` or raises `UncalibratedDataError`. |
| P5 | **Provenance ships with every plot.** A screenshot must be reproducible from the metadata it shows. | Each `Figure` carries a `metadata` dict with file SHAs, calibration version, function names, and timestamp. |
| P6 | **Citations are non-negotiable.** Every formula, threshold, or default has a primary-source reference. | Enforced by `documentation-standards.mdc`; CI fails if a new physics function has no `References` block. |
| P7 | **Tests precede expansion.** A new feature does not land without a regression test. | The CI gate is hard; "test will come later" PRs are rejected. |
| P8 | **Streamlit is presentation, not logic.** No physics or I/O lives in `app.py`. | `app.py` may only call into `blackhole.*` modules and render. Any computation inside `app.py` is a smell. |

---

## 2. Architecture additions

These are the new modules and contracts Phase 2 introduces. They are
declarative; each is implemented in one or more milestones (see §3 and §4).

### 2.1 Source catalog (`blackhole/catalog.py`)

```python
@dataclass(frozen=True)
class Source:
    name: str                    # canonical name, e.g. "NGC 1068"
    aliases: tuple[str, ...]     # SIMBAD-resolvable alternatives
    coord: SkyCoord              # frame=ICRS
    redshift: float | None
    distance_mpc: float | None   # luminosity distance
    distance_ref: str            # "Tully+2013", "Riess+2022", etc.
    m_bh_msun: float | None
    m_bh_err_msun: float | None
    m_bh_ref: str
    type: Literal[
        "seyfert1", "seyfert2", "llagn", "xrb_hmxb", "xrb_lmxb",
        "lrd", "quasar", "blazar", "tde",
    ]
    notes: str = ""
```

Catalog is a `tuple[Source, ...]` constant. **No mutation at runtime.**
Loaded by `from blackhole.catalog import CATALOG, by_short_id`.

### 2.2 Calibration registry (`blackhole/calibration.py`)

Each survey gets a callable
`(image: ImageData) -> CalibratedImage` where
`CalibratedImage.array` carries `astropy.units.Quantity` per pixel
(Jy/pix, MJy/sr, erg/s/cm²/keV, etc.) plus the original `ImageData`.

```python
@dataclass(frozen=True)
class CalibratedImage:
    array: u.Quantity            # 2D, with units
    pixel_solid_angle: u.Quantity  # steradian/pixel
    zeropoint: float | None
    zeropoint_ref: str
    method: str                  # "MAGZP+EXPTIME", "BUNIT", "ECF", ...
    original: ImageData

class UncalibratedDataError(RuntimeError): ...

CALIBRATORS: dict[str, Callable[[ImageData], CalibratedImage]] = {
    "2MASS-K": _calibrate_2mass,
    "DSS":     _calibrate_dss,    # raises UncalibratedDataError
    "IRIS_12": _calibrate_iris,
    ...
}

def calibrate(image: ImageData) -> CalibratedImage:
    """Dispatch on header (TELESCOP, INSTRUME, SURVEY, BUNIT)."""
```

DSS, plate-scan data, and any source without a documented calibration
**raises**, never silently returns a number.

### 2.3 Photometry contract (`blackhole/photometry.py`)

```python
def aperture_photometry(
    cal_image: CalibratedImage,
    source: Source,
    aperture_radius: u.Quantity,            # arcsec
    background_annulus: tuple[u.Quantity, u.Quantity],  # (r_in, r_out)
    *,
    method: Literal["exact", "subpixel"] = "exact",
) -> SEDPoint:
    """Returns SEDPoint with .flux_density and .flux_err set, units consistent."""
```

Implemented on top of `photutils`. Returns an `SEDPoint` with `band`,
`source`, `flux_err`, and `upper_limit` (set to True if measured flux <
3·σ_background).

### 2.4 Cosmology service (`blackhole/cosmo.py`)

Thin wrapper around `astropy.cosmology.Planck18` exposing:

```python
def luminosity_distance(z: float) -> u.Quantity
def angular_diameter_distance(z: float) -> u.Quantity
def kpc_per_arcsec(z: float) -> u.Quantity
def k_correction(rest_lambda: u.Quantity, obs_lambda: u.Quantity,
                 z: float, alpha: float = -0.5) -> float
```

Hardcoded to Planck18 with a `set_cosmology(...)` escape hatch. Justification
inline in module docstring.

### 2.5 Provenance & version manifest (`blackhole/provenance.py`)

```python
@dataclass(frozen=True)
class Provenance:
    fits_sha256: str
    fits_path: str
    calibration_version: str    # CALIBRATORS table version
    function_chain: tuple[str, ...]  # ("calibrate", "aperture_photometry")
    library_version: str        # blackhole.__version__
    timestamp_utc: str
```

`Figure` instances get `.metadata["provenance"] = Provenance(...)`.
The UI's *Provenance* expander reads from this. Every saved figure also gets
a sidecar JSON.

### 2.6 Cache layer

- `@st.cache_resource` for `load_events`, `load_image`, `calibrate`.
- `joblib.Memory(location=".joblib_cache")` for `aperture_photometry` keyed
  on `(fits_sha256, source_name, aperture_arcsec)`.
- Cache invalidation: `BLACKHOLE_VERSION` env var bump, or `--force` in CLI.

---

## 3. Workstreams and milestone dependency graph

Eleven milestones, three pillars. Pillar A unblocks everything; Pillar B is
Phase 2 proper; Pillar C is the citable-artifact bar.

```
PILLAR A — Close the loop
─────────────────────────
M0  Regression net & repo hardening   (no deps)
M1  Source catalog                    (M0)
M2  Calibration registry              (M0)
M3  Aperture photometry pipeline      (M1, M2)
M4  Light-curve correctness (GTI)     (M0)
M5  Spectrum tab honesty              (M0; opt. Sherpa hook = M5b)
M6  Provenance panel & UI revamp      (M3, M4, M5)

PILLAR B — Phase 2: Little Red Dots
───────────────────────────────────
M7  LRD catalog & cosmology           (M1)
M8  JWST ingestion                    (M2)
M9  LRD diagnostics                   (M3, M7, M8)
M10 Composite imagery (RGB)           (M2)

PILLAR C — Citable artifact
───────────────────────────
M11 Worked notebooks + Zenodo release (M3, M6, M9)
```

A milestone is "done" only when its **exit criteria** all pass.

---

## 4. Per-milestone detail

Each milestone uses the same template. Treat exit criteria as binary; treat
"out of scope" as a contract.

---

### UX iteration 1 (between M0 and M1) — Overview tab, in-page selector, per-tab controls

**Status.** Landed on `main` (branch: `ux-overview-and-controls`).

**Scope.**
- New **Overview** tab — thumbnail grid of every FITS in `fits_data/`,
  rendered through the same `wcs_plot` / `bin_to_image` paths the rest
  of the UI uses; cached on (path, stretch, cmap).
- File picker moved from the sidebar into a top-of-page selector
  alongside the metric strip.
- Stretch + cmap moved from the sidebar into per-tab **Display options**
  expanders on the Image and X-ray Events tabs.
- SED tab auto-detects the target from the active filename
  (`ngc1068`, `m87`, `cygx1` patterns); manual override remains.
- Sidebar trimmed to orientation content only.

**Deferred to M1 (catalog-driven UX).**
- Target-aware "science banner" above the fold (Eddington luminosity at
  the catalog mass, expected radio/IR/X-ray loudness ratios, classification
  call-outs). Cannot ship cleanly until M1 makes `BHSource` a real
  object with mass, distance, classification.

### M0 — Regression net & repo hardening

**Goal.** Build the safety net that lets every subsequent milestone be merged
with confidence. Nothing else lands until this does.

**Scope.**
- Convert layout to `pyproject.toml` (project metadata, `[project.scripts]`
  entry for `blackhole-download = scripts.download_data:main`).
- Add `tests/` with `pytest` + `pytest-cov`.
- Add `.github/workflows/ci.yml`: matrix on Python 3.11, 3.12; jobs:
  `ruff check`, `ruff format --check`, `mypy --strict blackhole/`,
  `pytest -q --cov=blackhole --cov-fail-under=70`.
- Add `ruff.toml` and `mypy.ini`.
- Move loose `.DS_Store` and `__pycache__` out of git; harden `.gitignore`.
- Add `pip-compile`-generated `requirements-lock.txt` for byte-reproducible
  installs.
- Create `tests/conftest.py` with a `tiny_fits` fixture (synthesized 64×64
  image with WCS, 1000-row event-list HDU, GTI extension).
- Seed tests for every existing module:
  - `test_io.py`: `inspect()` finds correct HDU count; `load_image()` returns
    proper shape and WCS-detection; `load_events()` parses our synthetic
    EVENTS HDU; `bin_to_image()` produces non-zero counts for sample events.
  - `test_sed.py`: Jy ↔ erg/s/cm²/Hz ↔ νFν round-trip across the four band
    helpers; `to_nu_fnu()` matches a hand-computed value at 1.4 GHz, 1 Jy
    to within 1e-12.
  - `test_lightcurves.py`: Lomb-Scargle on a 100-Hz sine recovers the
    frequency to within Δf.
  - `test_physics.py`:
    - `eddington_luminosity(1.0) == 1.26e38 erg/s` (4 sig figs).
    - `eddington_ratio(L_Edd, 1.0) == 1.0` exactly.
    - `accretion_rate_msun_yr(1e44, 0.1)` matches hand value.
    - `hardness_ratio(0, 100) == 1.0`; `hardness_ratio(100, 0) == -1.0`.
    - `donley_2012_agn` corners match published Fig. 1 region.
    - `fractional_rms` on noise-only data returns NaN.
    - `fractional_rms` on a known sinusoid recovers σ/⟨x⟩ within 5%.

**Exit criteria.**
- [x] `pytest -q --cov` reports ≥70% coverage across `blackhole/`.
      **Achieved: 92% (100 tests passing).**
- [x] `ruff check` passes.
- [x] `mypy --strict blackhole/` passes.
- [x] GitHub Actions CI green on `main` (workflow merged; runs in PR).
- [x] `pip install -e ".[dev]"` works.
- [x] No `.DS_Store` or `__pycache__` in `git status`.

**Files created.** `pyproject.toml`, `ruff.toml`, `mypy.ini`,
`.github/workflows/ci.yml`, `tests/conftest.py`, `tests/test_io.py`,
`tests/test_sed.py`, `tests/test_spectra.py`, `tests/test_lightcurves.py`,
`tests/test_wcs_plot.py`, `tests/test_physics/test_accretion.py`,
`tests/test_physics/test_infrared.py`, `tests/test_physics/test_variability.py`,
`tests/test_physics/test_spectral_xray.py`.

**Deferred from M0 (tracked, not blocking).**
- `ruff format --check` in CI. The first run produces ~1740 lines of
  cosmetic churn across 18 files; mixing that with M0's substantive
  changes would obscure the test/CI scaffolding. **Action**: a separate
  follow-up PR (`chore/ruff-format-baseline`) will run `ruff format .`
  and add `ruff format --check` to the CI workflow.
- `pip-compile`-generated `requirements-lock.txt`. Adds operational
  surface (lockfile regeneration cadence, dependabot config) that we
  haven't designed yet. **Action**: introduced with M5 (caching &
  reproducibility) when we already touch determinism plumbing.

**Bug fixed during M0.**
- `blackhole.io.load_image`: previously returned a non-None WCS for FITS
  files with no celestial CTYPE keywords (astropy quietly returns an
  identity 2-axis WCS). Now requires `WCS.has_celestial`. Matches
  `docs/PITFALLS.md` #3 and `04-fits-handling.mdc`.

**Risks.**
- Synthetic FITS in `conftest.py` drifting from reality: TODO in M1 — diff
  fixture-header keys against one real cutout in `fits_data/` after the
  source-catalog work makes real cutouts trivially loadable.

**References.** astropy testing guide; pytest-cov docs; ruff config.

---

### M1 — Source catalog

**Goal.** A single typed catalog of every studied source, with physical
parameters that the rest of the codebase reads instead of re-deriving.

**Scope.**
- Implement `blackhole/catalog.py` per §2.1.
- Seed entries for NGC 1068, M87, Cyg X-1 with peer-reviewed values:

  | Field | NGC 1068 | M87 | Cyg X-1 |
  |---|---|---|---|
  | M_BH (M☉) | (1.66 ± 0.04)×10⁷ | (6.5 ± 0.7)×10⁹ | 21.2 ± 2.2 |
  | M_BH ref | Lodato+2003 (maser) | EHT Coll. 2019 Paper VI | Miller-Jones+2021 |
  | z / d | z=0.00379; d=14.4 Mpc | d=16.8 Mpc (Bird+2010) | d=2.22 kpc |
  | Type | seyfert2 | llagn | xrb_hmxb |

- Helper functions: `by_short_id`, `resolve_simbad`, `eddington_luminosity_of`,
  `distance_to`, `redshift_of`.
- UI metadata banner reads from the catalog: shows M_BH, d, type, L_Edd at
  the top of every selected target's view.

**Exit criteria.**
- [ ] `tests/test_catalog.py` checks every numerical field against its
      reference paper value (≤1% drift) for all three Phase-1 targets.
- [ ] `eddington_luminosity_of("ngc1068")` matches a hand-computed
      L_Edd ≈ 2.1×10⁴⁵ erg/s.
- [ ] `app.py` imports `from blackhole.catalog import ...`; selecting a file
      whose filename prefix matches a `short_id` displays its catalog entry
      in the banner.
- [ ] No magic numbers (M_BH, distance, redshift) remain in `app.py` or
      `blackhole/sed.py`.

**Files created.** `blackhole/catalog.py`, `tests/test_catalog.py`.
**Files modified.** `app.py` (banner section).

**Risks.** Some published M_BH values disagree by factors of 2–10
(NGC 1068's pre-2003 estimates). Pin to a single primary source per entry
and document the disagreement in `notes`.

**References.** SIMBAD; NED; Lodato & Bertin 2003; EHT Collaboration 2019,
ApJL 875, L6; Miller-Jones+2021, Science 371, 1046.

---

### M2 — Calibration registry

**Goal.** Make every supported survey return calibrated, units-bearing data;
make every unsupported survey raise loudly.

**Scope.**
- Implement `blackhole/calibration.py` per §2.2.
- Implement calibrators for the surveys in our current dataset:

  | Survey | Method |
  |---|---|
  | **2MASS-K** | flux = 10^((MAGZP - mag)/2.5) · F_ν,0(K); F_ν,0 = 666.7 Jy (Cohen+2003) |
  | **IRIS 12/25/60/100** | BUNIT = MJy/sr; multiply by pixel solid angle |
  | **AKARI WIDE-S (90 µm)** | BUNIT = MJy/sr |
  | **VLA FIRST 1.4 GHz** | BUNIT = Jy/beam; deconvolve via BMAJ/BMIN |
  | **RASS broad** | EXPTIME × ECF (0.1–2.4 keV) from Snowden+1997 |
  | **DSS** | raises `UncalibratedDataError` |
- Each calibrator's docstring must include the primary-source reference,
  the BUNIT it expects, and what it returns. (Enforced by
  `documentation-standards.mdc`.)
- One unit test per calibrator using a hand-built synthetic FITS with known
  pixel values and a known answer.

**Exit criteria.**
- [ ] `tests/test_calibration.py`: for each survey, a 2-pixel synthetic image
      converts to the expected `Quantity` to within 1e-4.
- [ ] `calibrate(dss_image)` raises `UncalibratedDataError`.
- [ ] `mypy --strict` passes on the new module.
- [ ] Every public function has a `References` block in its docstring.

**Files created.** `blackhole/calibration.py`, `tests/test_calibration.py`,
`docs/calibration_table.md`.

**Risks.**
- IRIS pixel solid angle depends on CDELT1/CDELT2 — must read from header
  *per file*, not assume.
- VLA FIRST beam is asymmetric; for unresolved sources, peak-flux-density
  in Jy/beam *is* the flux. For extended emission this is a lower limit.
  Document this trap explicitly.

**References.** Cohen, Wheaton & Megeath 2003 AJ 126, 1090; Wright+2010
(WISE); Miville-Deschênes & Lagache 2005 (IRIS); Becker+1995 (FIRST);
Snowden+1997 (RASS ECF).

---

### M3 — Aperture photometry pipeline

**Goal.** The SED tab's numbers come from the local FITS files. The hardcoded
literature dict shrinks to a comparison overlay, not the data source.

**Scope.**
- Add `photutils` to requirements.
- Implement `blackhole/photometry.py` per §2.3.
- Background estimation: median-clipped annulus, σ from MAD.
- Upper-limit detection: if measured flux < 3σ_bg, return `SEDPoint` with
  `upper_limit=True` and `flux_density = 3σ_bg`.
- Update SED tab: build `SED` from photometry calls, with a checkbox
  "Overlay published literature values" that draws the old hardcoded table
  as faint comparison markers.
- Add aperture-size sliders per band group in the UI (radio uses much larger
  apertures than optical).

**Exit criteria.**
- [ ] `tests/test_photometry.py`: synthetic point source with known
      total flux → recovered flux within 2% (when PSF-fit) or within
      aperture-correction tolerance (when raw).
- [ ] For NGC 1068, the 2MASS K-band aperture flux from the local cutout
      matches the 2MASS XSC catalog value to within 10% (smoke test).
- [ ] No literal flux values remain inside `app.py`. The literature dict
      is moved to `docs/literature_seds.py` and is only used for overlay.
- [ ] UI shows aperture radius, background annulus, and 1-σ uncertainty for
      every plotted point.

**Files created.** `blackhole/photometry.py`,
`tests/test_photometry.py`, `docs/literature_seds.py`.
**Files modified.** `app.py` (SED tab), `requirements.txt` (add `photutils`).

**Risks.**
- M87's extended halo and jet make "aperture photometry" ill-defined.
  Either expose multiple aperture choices or restrict M87 to a core
  measurement and document this prominently.
- NGC 1068 has bright satellite knots in the optical at small radii. The
  default DSS aperture may need a higher inner annulus.

**References.** Bradley+ (`photutils` docs); 2MASS Extended Source Catalog;
Jarrett+2003.

---

### M4 — Light-curve correctness (GTI)

**Goal.** A light curve from a real X-ray event file has no fake zero-count
gaps; F_var is computed only over valid exposure.

**Scope.**
- `blackhole/io.py`: extend `EventList` with `gti: np.ndarray | None`
  (shape (N, 2), columns `(start, stop)`).
- `load_events()` reads the GTI extension (HDU named `GTI` or `STDGTI` or
  per-CCD GTI) into that field.
- `lightcurves.bin_events_to_lightcurve` accepts `apply_gti: bool = True`;
  bins overlapping no GTI are dropped from the output rather than rendered
  as zero.
- Effective exposure per bin is computed as the intersection of the bin and
  the GTIs; Poisson error uses effective exposure, not bin width.
- F_var uses only bins with `effective_exposure / bin_size > 0.5`.

**Exit criteria.**
- [ ] `tests/test_lightcurves.py::test_gti_masks_gaps`: synthetic event list
      with a 1 ks gap shows the bin during the gap is absent from output,
      not zero.
- [ ] F_var on a known constant source through a gapped GTI is consistent
      with zero (NaN, by the existing convention).
- [ ] The light-curve plot title shows the GTI exposure fraction.

**Files modified.** `blackhole/io.py`, `blackhole/lightcurves.py`,
`tests/test_lightcurves.py`.

**Risks.**
- Some missions split GTIs per detector/CCD; union them carefully.
- NuSTAR has per-FPMA/FPMB GTIs that are not always identical.

**References.** Chandra POG §4; XMM-Newton SAS users' guide GTI section;
Vaughan+2003 (F_var with exposure weighting).

---

### M5 — Spectrum tab honesty

**Goal.** Either remove the channel-space "power-law fit" or replace it with
a calibrated forward-folded fit. No middle ground.

**Scope (M5a — minimum).**
- Rename the displayed parameter from "Γ (photon index)" to
  "α_channel (descriptive)" and clearly label the axis "Channel (not energy)."
- Remove the `classify_photon_index(g)` call from the channel-space fit
  output. That label only applies to a calibrated Γ.

**Scope (M5b — Sherpa hook, optional but recommended).**
- Add `sherpa` to optional extras: `pip install BlackHoleResearch[xray]`.
- Implement `blackhole/sherpa_bridge.py` that, given a PHA + RMF + ARF
  triplet, runs a `xs_powerlaw + xs_phabs` fit and returns
  `(Gamma, Gamma_err, NH, NH_err, chi2_reduced)`.
- Add a "Calibrated fit (Sherpa)" toggle in the Spectrum tab. When response
  files are present in `fits_data/`, the channel-space curve disappears and
  the Sherpa fit shows instead.
- Ship one example Chandra ACIS-S spectrum (`ngc1068_acis.pha`) with its
  RMF/ARF for the demo.

**Exit criteria (M5a).**
- [ ] No label in the UI claims a photon-index value when only channel-space
      data is available.
- [ ] `tests/test_spectra.py::test_no_photon_index_without_rmf` asserts
      that `classify_photon_index` is *not* called from the channel-space
      code path.

**Exit criteria (M5b, if pursued).**
- [ ] `tests/test_sherpa_bridge.py`: fitting a Sherpa-generated fake spectrum
      with input Γ=2.0 recovers Γ within fit uncertainties.
- [ ] Documentation in `docs/spectral_fitting.md` describes how to install
      Sherpa and what RMF/ARF to use for each shipped example.

**Files modified.** `app.py` (Spectrum tab), `blackhole/spectra.py`
(remove/rename `classify_*` call), possibly `blackhole/sherpa_bridge.py`.

**Risks.**
- Sherpa's installation footprint is non-trivial. Keep as an *optional* extra
  so the core install stays slim.
- Forward-folding fits are slow; cache by SHA.

**References.** OGIP/92-007; CIAO Sherpa documentation; Freeman+2001.

---

### M6 — Provenance panel & UI revamp

**Goal.** Every plot in the UI is reproducible from on-screen metadata.

**Scope.**
- Implement `blackhole/provenance.py` per §2.5.
- All `render_*` functions attach a `Provenance` to `fig.metadata`.
- Add a "Provenance" expander under each tab that renders the metadata as a
  table, plus a "Copy citation" button that emits the right BibTeX snippet.
- Refactor dark-mode styling into one `blackhole/_style.py::apply_dark(ax)`
  helper, used by every renderer.
- Add a sidebar checkbox "Save figure + sidecar JSON to `outputs/`".
- Banner shows: target, M_BH, d, z, L_Edd, type.

**Exit criteria.**
- [ ] Every tab in the app has a Provenance expander.
- [ ] `tests/test_provenance.py`: saving a figure produces a `.png` and a
      `.json` sidecar; loading the JSON reproduces the `Provenance` struct.
- [ ] No `for spine in ax.spines.values(): spine.set_color("white")` appears
      more than once in the codebase (grep test).

**Files created.** `blackhole/provenance.py`, `blackhole/_style.py`,
`tests/test_provenance.py`.
**Files modified.** every `blackhole/*.py` with a `render_*` function.

**Risks.** None significant. Largely refactoring.

---

### M7 — LRD catalog & cosmology

**Goal.** A Phase-2 set of LRD targets in the catalog with cosmologically
correct distances, ready for diagnostics.

**Scope.**
- Implement `blackhole/cosmo.py` per §2.4.
- Add LRD entries to the catalog from published samples:

  | Source | z | sample | reference |
  |---|---|---|---|
  | CEERS-415 | 5.62 | CEERS | Harikane+2023 |
  | CEERS-2782 | 5.24 | CEERS | Harikane+2023 |
  | JADES-GN-954880 | 5.49 | JADES | Maiolino+2024 |
  | COS-XL-747 | 7.04 | COSMOS-Web | Akins+2025 |
  | RUBIES-EGS-49140 | 6.68 | RUBIES | Wang+2024 |

  For each: redshift, broad-Hα FWHM (km/s), measured M_BH (from Reines &
  Volonteri 2015 or Greene & Ho 2005 calibration), L_X 2 keV upper limit,
  F_2500 if measured, references.
- `type = "lrd"`.

**Exit criteria.**
- [ ] `tests/test_cosmo.py`: `luminosity_distance(1.0)` returns
      6701 Mpc (Planck18) ± 0.5 Mpc.
- [ ] `tests/test_catalog.py` validates every LRD entry has a non-None
      `redshift` and a `distance_mpc` matching Planck18 to ≤0.1%.
- [ ] At least 5 LRDs in the catalog with primary-source citations.

**Files created.** `blackhole/cosmo.py`, `tests/test_cosmo.py`.
**Files modified.** `blackhole/catalog.py`.

**References.** Planck Collaboration 2020 (Planck18); Greene+2024;
Maiolino+2024; Harikane+2023; Akins+2025; Wang+2024.

---

### M8 — JWST ingestion

**Goal.** Read JWST NIRSpec 1D spectra and NIRCam photometry into the same
data model used for the Phase-1 cutouts.

**Scope.**
- `blackhole/jwst.py` handles:
  - **NIRSpec PRISM/G395M `_x1d.fits`**: returns `specutils.Spectrum1D` with
    wavelength in microns, flux in Jy, uncertainty in Jy.
  - **NIRCam `_i2d.fits`** (calibrated mosaic): returns `CalibratedImage`
    with proper photometric calibration via PHOTMJSR keyword.
  - **NIRCam `_cat.fits`** source catalogs: returns an astropy `Table`.
- Add a NIRSpec viewer tab to the Streamlit UI: line markers for Hα, Hβ,
  [O III] 5007, [O III] 4959, [N II] 6583, He II 4686, [S II] 6716/6731,
  Mg II 2798, C IV 1549, Lyα 1216, all redshift-aware.
- Add a `download_jwst.py` extension to `scripts/` using `astroquery.mast`
  for the LRDs in the catalog.

**Exit criteria.**
- [ ] `tests/test_jwst.py`: synthetic NIRSpec-like `_x1d.fits` loads and
      returns a `Spectrum1D` of correct length and units.
- [ ] For one LRD with shipped public NIRSpec data, the redshift fit on
      [O III] 5007 matches the published z to within 0.01.
- [ ] NIRSpec viewer renders rest-frame line markers correctly at z=5 in
      a unit test using `pytest-mpl` baseline images.

**Files created.** `blackhole/jwst.py`, `scripts/download_jwst.py`,
`tests/test_jwst.py`.
**Files modified.** `app.py` (new tab), `requirements.txt` (add
`specutils`, already present).

**Risks.**
- MAST queries are slow and rate-limited. Cache aggressively.
- JWST data products evolve between pipeline versions; pin to a specific
  CRDS context in `download_jwst.py` and document.

**References.** Böker+2023 (NIRSpec); Rieke+2023 (NIRCam);
JWST Data Handbook §3.1 and §4.4.

---

### M9 — LRD diagnostics

**Goal.** Reproduce the three canonical LRD diagnostic plots from local data.

**Scope.**
- **D1 — L_X / L_UV vs Lusso & Risaliti 2016.** Compute monochromatic L at
  rest-frame 2 keV (from X-ray data, or upper limit) and 2500 Å (from
  NIRSpec). Plot against the LR16 mean relation
  `log L_2keV = 0.642 log L_2500 + 7.07`. Mark LRDs with α_ox < expected by
  ≥0.3 (they all are).
- **D2 — Balmer break.** Fit a piecewise power-law `f_λ ∝ λ^α₁` for
  3500–3645 Å rest and `f_λ ∝ λ^α₂` for 3700–4500 Å. Report break amplitude
  in magnitudes. Compare to Setton+2024 / Wang+2024 LRD distribution.
- **D3 — Broad-Hα + narrow [O III] fit → M_BH.** Gaussian decomposition
  (broad Hα + narrow Hα + [N II] doublet). Apply Reines & Volonteri 2015:
  `log M_BH = log ε + 6.57 + 0.47 log(L_Hα / 1e42) + 2.06 log(FWHM_Hα / 1000)`.
  Compare to Greene+2024 published M_BH for the same source.
- **D4 — Donley wedge with LRDs overlaid.** Use the existing `donley_2012_agn`;
  scatter the LRDs in the W1−W2 vs W2−W3 plane (with upper-limit arrows
  for non-detections, which is most of them).

Each diagnostic gets a tab or sub-tab.

**Exit criteria.**
- [ ] `tests/test_lrd_diagnostics.py`:
  - LR16 mean relation matches Eq. (3) of LR16 at two test points.
  - Balmer break detector recovers a synthetic 0.5-mag break to within
    0.05 mag.
  - Reines-Volonteri recovers M_BH for a synthetic broad line of known
    L_Hα and FWHM to within 0.1 dex.
- [ ] For one published LRD with reported M_BH (e.g. JADES-GN-954880),
      the pipeline value matches Maiolino+2024 within their stated 0.3 dex
      systematic uncertainty.
- [ ] All four diagnostics render in the Streamlit UI with provenance
      panels.

**Files created.** `blackhole/lrd.py`, `tests/test_lrd_diagnostics.py`.
**Files modified.** `app.py`.

**Risks.**
- Reines & Volonteri 2015 calibration was derived for low-z dwarfs; its
  applicability to z>5 LRDs is itself a research question. Surface this
  in the UI as a footnote citing Maiolino+2024 §5 on the calibration
  question.

**References.** Lusso & Risaliti 2016 ApJ 819, 154 Eq. (3);
Reines & Volonteri 2015 ApJ 813, 82; Setton+2024; Wang+2024 ApJ 969, 13;
Greene+2024 ApJ 964, 39; Maiolino+2024.

---

### M10 — Composite imagery (RGB) — optional polish

**Goal.** Lupton+2004 asinh RGB composites for any three-band combination.

**Scope.**
- `blackhole/_rgb.py::lupton_rgb(r, g, b, *, Q, stretch)` based on
  `astropy.visualization.make_lupton_rgb` with re-projection via
  `reproject` to align inputs that have different WCS pixel scales.
- A "Composite" tab in the UI allowing the user to choose R/G/B band
  assignments from any three loaded files for the current target.

**Exit criteria.**
- [ ] `tests/test_rgb.py`: a flat RGB of three constant images yields a
      uniform gray output.
- [ ] M87 RGB (optical/IR/radio) renders without WCS-alignment errors.

**Files created.** `blackhole/_rgb.py`, `tests/test_rgb.py`.
**Files modified.** `app.py`, `requirements.txt` (add `reproject`).

**References.** Lupton, Blanton, Fekete+ 2004 PASP 116, 133.

---

### M11 — Worked notebooks + Zenodo release

**Goal.** The DOI bar from §0.5 — turn this into something that can be cited
in a paper.

**Scope.**
- `docs/notebooks/01_ngc1068_seyfert2.ipynb`: end-to-end from `download_data`
  to SED with calibrated 2MASS-K + IRIS aperture photometry; reproduce the
  characteristic SED with the mid-IR torus bump; show consistency with
  published values.
- `docs/notebooks/02_m87_jet.ipynb`: VLA + DSS + ROSAT composite; aperture
  photometry of the core; Eddington-ratio calculation.
- `docs/notebooks/03_cygx1_xray.ipynb`: light-curve + F_var + GTI handling
  on a Chandra event file.
- `docs/notebooks/04_lrd_diagnostics.ipynb`: full Phase-2 pipeline on one
  published LRD; reproduces published M_BH and α_ox within stated tolerances.
- Tag `v0.2.0`; archive on Zenodo; add DOI badge to README.
- Create `CITATION.cff` and `CITATION.md`.

**Exit criteria.**
- [ ] All four notebooks execute top-to-bottom on CI (`nbmake` plugin) in
      under 5 minutes total.
- [ ] DOI present in `README.md` and `CITATION.cff`.
- [ ] Zenodo upload includes the full repo at the tagged commit.

**Files created.** `docs/notebooks/*.ipynb`, `CITATION.cff`, `CITATION.md`.
**Files modified.** `README.md`, `.github/workflows/ci.yml` (add `nbmake`
job).

**References.** Zenodo–GitHub integration; Citation File Format spec.

---

## 5. Phase 2 deep dive — Little Red Dots

This section is the *scientific* content of Phase 2. It is normative for
M7–M9; everything below must be cited in the relevant module docstring.

### 5.1 Background

JWST imaging in 2023–2024 (CEERS, JADES, RUBIES, COSMOS-Web) revealed a
population of point-like, compact, red-continuum sources at z ≈ 4–9 with:

- **Broad permitted lines** (broad Hα FWHM ≈ 1000–4000 km/s) implying gas
  velocities consistent with accretion onto a 10⁶–10⁸ M☉ black hole. This
  is the unambiguous AGN signature.
- **No detected X-ray emission**, even in deep Chandra stacks. Akins+2025
  and Yue+2024 report 2σ upper limits ≳ 2 dex below the expected L_X for
  an AGN of inferred L_bol.
- **No detected hot-dust mid-IR torus**, even in MIRI stacks
  (Akins+2025 §4).
- **A characteristic v-shaped UV+optical SED** with a blue UV continuum
  and a red optical continuum, often interpreted as a Balmer break
  (Setton+2024, Wang+2024) indicating an evolved stellar population, or as
  reddened AGN continuum (debated).

### 5.2 Why this matters

The AGN unification picture (Antonucci 1993; Urry & Padovani 1995) requires
**every** accreting BH to power a hot corona (X-ray) and to be surrounded
by a dusty torus (mid-IR). LRDs appear to violate both. Either:

1. They are super-Eddington accretors with radiation-driven outflows that
   suppress the corona and disrupt the torus (Pacucci & Narayan 2025;
   Lambrides+2024 "buried AGN").
2. They are not classical AGN — possibly direct-collapse seed remnants or
   accreting BHs in unusual host environments (Maiolino+2024; Inayoshi+2025).
3. They are heavily Compton-thick along *all* lines of sight, with column
   densities ≳ 10²⁵ cm⁻² (Greene+2024 §5).

The diagnostics in M9 are designed to test these hypotheses.

### 5.3 Target list (seed)

See §M7 table. Add as JWST cycles continue.

### 5.4 Diagnostics — exact definitions

- **α_ox = −0.3838 · log(F_2keV / F_2500Å)** (Tananbaum+1979).
- **L_2keV expected from L_2500Å**: LR16 Eq. (3),
  `log L_2keV = 0.642 log L_2500 + 7.07`, scatter 0.24 dex.
- **Balmer break amplitude**: 2.5 · log(⟨f_λ(3700–4500)⟩ / ⟨f_λ(3300–3700)⟩),
  computed in rest frame.
- **M_BH (Hα)**: Reines & Volonteri 2015 Eq. (5),
  `log M_BH = 6.57 + 0.47 log(L_Hα / 1e42) + 2.06 log(FWHM_Hα / 1000)`.

### 5.5 Open research questions to flag in UI

- Is the inferred M_BH from broad Hα trustworthy when the standard
  reverberation-mapping calibration was derived at z < 0.3?
- Are the X-ray non-detections real Compton-thick obscuration or a sign of
  a fundamentally different accretion state?
- Is the Balmer break stellar or AGN?

These should appear as expandable "Open question" boxes alongside the
diagnostic plots.

---

## 6. Hardening (cross-cutting; revisited in every PR)

### 6.1 Test taxonomy
- **Unit tests.** One per public function. Synthetic inputs, hand-computed
  outputs.
- **Regression tests.** "Physics constants" — Eddington luminosity, Stefan-
  Boltzmann, etc. — must reproduce textbook values to specified precision.
- **Integration tests.** End-to-end: download (mock) → calibrate → photometry
  → SEDPoint → render.
- **Notebook tests.** `pytest --nbmake docs/notebooks/`.
- **Image regression tests (optional).** `pytest-mpl` baseline images for
  every renderer.

### 6.2 CI matrix
- Python: 3.11, 3.12.
- OS: ubuntu-latest, macos-latest.
- Jobs: `ruff check`, `ruff format --check`, `mypy --strict blackhole/`,
  `pytest --cov`, `nbmake docs/notebooks/`, `pip install -e .` dry-run.

### 6.3 Coverage targets
- Overall: ≥80%.
- `blackhole/photometry.py`, `blackhole/calibration.py`,
  `blackhole/lrd.py`: ≥90%.
- `app.py`: not measured (UI).

### 6.4 Type discipline
- Public functions: full annotations, including return types.
- Generic numerics: prefer `np.ndarray` and `astropy.units.Quantity` over
  `Any`.
- `mypy --strict` runs on `blackhole/` (not `scripts/` or `app.py`).

### 6.5 Documentation discipline
Enforced by `.cursor/rules/documentation-standards.mdc`. Every module-level
docstring must include:
1. **What this module exists for** (one sentence).
2. **Key physics or convention** (3–10 lines).
3. **References** (primary sources, with DOI or arXiv ID where possible).
4. **Pitfalls** if the module touches one of the items in `PITFALLS.md`.

Every public function docstring must include:
- One-line summary.
- Parameter and return types (in addition to annotations).
- A `References` block when the function implements a published formula.

### 6.6 Doc build
- `mkdocs` with `mkdocs-material` theme.
- Auto-generated API reference via `mkdocstrings[python]`.
- `mkdocs gh-deploy` on every tag.

### 6.7 Release process
1. `pytest -q` green, coverage ≥80%.
2. `CHANGELOG.md` updated under "Unreleased."
3. `__version__` bumped in `blackhole/__init__.py`.
4. Tag `vX.Y.Z`; push; GitHub Release auto-creates from the tag.
5. Zenodo mints a DOI (GitHub–Zenodo integration).
6. README badge updated.

---

## 7. Phase 3 hints (record now, do later)

- **Sherpa / PyXspec full integration.** RMF/ARF handling, multi-component
  models, MCMC error analysis.
- **IXPE polarimetry.** Stokes Q/U handling.
- **Coordination with time-domain alerts** (e.g. ZTF, ATLAS, LSST broker).
- **EHT / VLBI ring fitting** for M87-class targets.
- **Galactic-extinction de-reddening** (Schlafly & Finkbeiner 2011) — small
  but real for many extragalactic sources.
- **K-corrections** beyond the simple power-law in §2.4.
- **Pulsar-timing array data** as a counterweight: not BH accretion, but
  shares the FITS/HDF5 infrastructure.

---

## 8. Decision log (ADR-style)

Each architectural decision is recorded in `docs/adr/`.

### ADR-0001 — Adopt typed Source catalog (M1)
**Status.** Accepted.
**Context.** Hardcoded numbers in `app.py` for M_BH and distance were the
biggest scientific-honesty risk. Different parts of the codebase used
inconsistent or absent values.
**Decision.** Centralize in `blackhole/catalog.py` as an immutable,
typed dataclass tuple.
**Consequences.** Adding a target is now a one-place change; consistency is
testable; tests can pin to primary-source values.

### ADR-0002 — Calibration is dispatch on header (M2)
**Status.** Accepted.
**Context.** Each survey has its own flux convention; doing this implicitly
in many places leads to drift.
**Decision.** Single `calibrate(image)` entry point that dispatches via a
registry keyed on `TELESCOP`/`INSTRUME`/`SURVEY`/`BUNIT`. Uncalibrated data
raises.
**Consequences.** Adding a survey is one new function + one registry entry +
one test. Misuse is a loud failure, not a silent wrong number.

### ADR-0003 — Streamlit is presentation, not logic (P8)
**Status.** Accepted.
**Context.** `app.py` was beginning to acquire physics (SED lookup tables).
**Decision.** Physics and I/O live in `blackhole/`; `app.py` may only call
into the library and render. A test grep-fails the build if `import numpy`
or `import astropy.units` appears in `app.py`.
**Consequences.** Refactoring the UI does not risk physics regressions;
the library is usable without Streamlit (notebooks, batch jobs).

### ADR-0004 — Channel-space "Γ" is removed or renamed (M5)
**Status.** Proposed (decided in M5a).
**Context.** Calling a channel-space fit's slope "Γ" sets up a false
comparison to the published photon index, which is energy-space and
forward-folded.
**Decision.** Rename to `α_channel` and remove the AGN-class label.
**Consequences.** Honesty preserved; users who want a real Γ are pointed
to the M5b Sherpa hook.

### ADR-0005 — Hard pin to Planck18 (M7)
**Status.** Proposed.
**Context.** Cosmology choice affects every L from D_L(z).
**Decision.** Use `astropy.cosmology.Planck18` by default; expose
`set_cosmology()` for advanced users; record the cosmology in `Provenance`.
**Consequences.** Reproducibility is guaranteed for the default; advanced
users have an escape hatch that is recorded in provenance.

---

## 9. Open research questions (researcher's wishlist)

These are *not* engineering tasks; they are scientific questions this tool
should make easier to investigate. They live here so future contributors can
see what would be a valuable application:

1. Are LRDs at fixed M_BH systematically less X-ray luminous than z~0 AGN
   matched in M_BH, or only at fixed L_bol? (Pipeline output: scatter plot.)
2. Do LRDs occupy a distinct region of the Donley wedge, or are they simply
   undetected? (Pipeline output: a stacked WISE/MIRI photometric upper-limit
   tally.)
3. Is the inferred Balmer break correlated with L_bol or with redshift?
   (Pipeline output: a `pandas.DataFrame` joinable with NIRSpec spectra.)
4. Among Phase-1 targets: how does M87's core flux at 1.4 GHz scale with
   aperture radius? (Pipeline output: aperture growth curve.)

---

## 10. Glossary

| Term | Definition |
|---|---|
| **AGN** | Active Galactic Nucleus — accreting supermassive black hole. |
| **ARF** | Auxiliary Response File — effective area vs energy. |
| **BUNIT** | FITS header keyword for the units of a pixel. |
| **EHT** | Event Horizon Telescope. |
| **F_var** | Fractional rms variability amplitude (Vaughan+2003). |
| **GTI** | Good Time Interval — when an instrument was accumulating valid data. |
| **HDU** | Header/Data Unit — one element of a FITS file. |
| **JWST** | James Webb Space Telescope. |
| **L_Edd** | Eddington luminosity, L = 4πGMm_p c / σ_T. |
| **LRD** | Little Red Dot — JWST-discovered red, compact, broad-line objects. |
| **NIRSpec** | Near-Infrared Spectrograph on JWST. |
| **OGIP** | Office of Guest Investigator Programs — the X-ray FITS conventions body. |
| **PHA** | Pulse Height Amplitude — X-ray spectrum FITS format. |
| **Planck18** | Planck Collaboration 2018 cosmology (Ω_m=0.315, H_0=67.4). |
| **RMF** | Response Matrix File — channel ↔ energy probability matrix. |
| **SED** | Spectral Energy Distribution — flux vs frequency on log-log axes. |
| **SkyView** | NASA cross-mission survey cutout service. |
| **WCS** | World Coordinate System — pixel ↔ sky-coordinate mapping. |
| **α_ox** | Optical-to-X-ray spectral index, −0.3838·log(F_2keV/F_2500Å). |
| **Γ** | Photon index for a power-law spectrum N(E) ∝ E^(−Γ). |

---

## 11. Appendix — primary-source references

Catalogued for citation in every module that uses them.

- **Antonucci, R., 1993, ARA&A 31, 473.** AGN unification.
- **Akins+2025, arXiv:2406.10341.** COSMOS-Web LRDs; IR/X-ray non-detections.
- **Becker+1995, ApJ 450, 559.** VLA FIRST survey description.
- **Bird+2010, MNRAS 405, 1409.** M87 distance.
- **Cohen, Wheaton & Megeath 2003, AJ 126, 1090.** 2MASS K-band F_ν,0.
- **Donley+2012, ApJ 748, 142.** Refined WISE AGN wedge.
- **Eddington, A. S., 1926.** *The Internal Constitution of the Stars.*
- **EHT Collaboration 2019, ApJL 875, L1–L6.** M87 event horizon.
- **Elvis+1994, ApJS 95, 1.** Mean quasar SED template.
- **Freeman+2001, ASP 238, 76.** Sherpa.
- **Greene & Ho 2005, ApJ 630, 122.** Hα-based BH mass.
- **Greene+2024, ApJ 964, 39.** LRD census in JADES.
- **Harikane+2023, ApJ 959, 39.** CEERS broad-line AGN.
- **Hickox & Alexander 2018, ARA&A 56, 625.** Obscured AGN review.
- **Inayoshi+2025.** LRDs as evolved direct-collapse remnants.
- **Lambrides+2024, ApJ 961, L25.** Buried-AGN interpretation of LRDs.
- **Lodato & Bertin 2003, A&A 398, 517.** NGC 1068 maser BH mass.
- **Lupton, Blanton, Fekete+ 2004, PASP 116, 133.** Asinh RGB.
- **Lusso & Risaliti 2016, ApJ 819, 154.** L_X / L_UV mean relation.
- **Maiolino+2024, A&A 691, A145.** JADES broad-line AGN catalog.
- **Mushotzky, Done & Pounds 1993, ARA&A 31, 717.** AGN X-ray review.
- **Miller-Jones+2021, Science 371, 1046.** Cyg X-1 BH mass refinement.
- **Miville-Deschênes & Lagache 2005, ApJS 157, 302.** IRIS reprocessing.
- **Novikov & Thorne 1973.** Thin disk GR derivation.
- **Pacucci & Narayan 2025.** Super-Eddington LRD model.
- **Planck Collaboration 2020, A&A 641, A6.** Planck18 cosmology.
- **Reines & Volonteri 2015, ApJ 813, 82.** Hα BH-mass calibration.
- **Ricci+2017, ApJS 233, 17.** BAT AGN spectroscopic survey.
- **Riess+2022, ApJ 934, L7.** Local H_0 distance ladder.
- **Setton+2024, ApJ 974, 4.** LRD Balmer-break interpretation.
- **Shakura & Sunyaev 1973, A&A 24, 337.** α-disk model.
- **Snowden+1997, ApJ 485, 125.** RASS energy conversion factors.
- **Stern+2012, ApJ 753, 30.** WISE single-color AGN cut.
- **Tananbaum+1979, ApJ 234, L9.** α_ox definition.
- **Urry & Padovani 1995, PASP 107, 803.** AGN unification synthesis.
- **VanderPlas 2018, ApJS 236, 16.** Lomb-Scargle astronomy primer.
- **Vaughan+2003, MNRAS 345, 1271.** F_var and excess variance.
- **van der Klis 2006.** Rapid X-ray variability review.
- **Wang+2024, ApJ 969, 13.** RUBIES LRDs; Balmer break analysis.
- **Webster & Murdin 1972, Nature 235, 37.** Cyg X-1 BH identification.
- **Wright+2010, AJ 140, 1868.** WISE Vega-AB conversions.
- **Yue+2024.** Deep Chandra stacks of LRDs.

---

*This document is normative for Phase 2. Changes to scope require an ADR
under §8. Changes to the citable-artifact criterion in §0.5 require a
maintainer decision recorded in the same place.*
