<div align="center">

# BlackHoleResearch

### A multiwavelength FITS explorer for black hole astrophysics — from ROSAT to JWST, all in one Streamlit app.

**Load NASA-archive FITS files from Chandra, XMM-Newton, NuSTAR, ROSAT, JWST (NIRSpec/NIRCam), WISE, 2MASS, VLA, and DSS. Render them correctly. Build calibrated multi-band SEDs. Bin X-ray event lists. Fit X-ray spectra. Compute Lomb–Scargle periodograms and Vaughan+2003 fractional rms variability. Reproduce Little Red Dot diagnostics. All on a laptop.**

[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg?logo=python&logoColor=white)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![astropy](https://img.shields.io/badge/powered%20by-astropy-orange.svg?logo=python)](https://www.astropy.org)
[![Streamlit](https://img.shields.io/badge/UI-Streamlit-FF4B4B.svg?logo=streamlit)](https://streamlit.io)
[![FITS](https://img.shields.io/badge/format-FITS-blueviolet.svg)](https://fits.gsfc.nasa.gov/)
[![NASA data](https://img.shields.io/badge/data-NASA%20public%20domain-005288.svg)](https://skyview.gsfc.nasa.gov/)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-D7FF64.svg)](https://docs.astral.sh/ruff/)
[![Status: Phase 1, Phase 2 planned](https://img.shields.io/badge/status-Phase%201%20released%20·%20Phase%202%20planned-success.svg)](PHASE2_PLAN.md)

[**Quick start**](#-quick-start)  ·
[**Features**](#-what-you-can-do-with-this)  ·
[**Targets**](#-the-targets-and-the-physics)  ·
[**Phase 2: Little Red Dots**](#-phase-2--jwst-little-red-dots-research-frontier)  ·
[**FITS pitfalls catalog**](docs/PITFALLS.md)  ·
[**Roadmap**](PHASE2_PLAN.md)  ·
[**Cite**](#-citation)

</div>

---

## What it does

A short list of the specific tasks this app is built around:

- Load a Chandra / XMM / NuSTAR event file and turn it into an image with the right energy cut. → [X-ray event list-to-image binning](#-x-ray-event-list-to-image-binning).
- View a JWST NIRSpec spectrum with redshifted line markers. → [Phase 2 NIRSpec tab (planned)](PHASE2_PLAN.md#m8--jwst-ingestion).
- Build a multi-band Spectral Energy Distribution for an AGN with consistent units across radio → X-ray. → [Unit-safe SED builder](#-multi-band-sed-builder-radio--mid-ir--optical--x-ray).
- Look up whether WISE colors are AB or Vega, and how the Stern / Donley cuts are defined. → [PITFALLS #6](docs/PITFALLS.md), [`physics/infrared.py`](blackhole/physics/infrared.py).
- Reproduce the L_X / L_UV plot for a Little Red Dot. → [Phase 2 diagnostics (planned)](PHASE2_PLAN.md#m9--lrd-diagnostics).
- Investigate zero-count gaps in an X-ray light curve. → [GTI handling (planned M4)](PHASE2_PLAN.md#m4--light-curve-correctness-gti).
- See the asinh stretch, the Donley wedge, the Eddington luminosity, the Lomb–Scargle periodogram, and the Vaughan+2003 F_var in one place with primary-source citations. → That's the whole app.

> **Search keywords (for discoverability):**
> Chandra event file Python viewer, JWST NIRSpec FITS Streamlit, multiwavelength SED AGN, OGIP PHA spectrum power law, WCS overlay matplotlib origin lower, asinh stretch astronomical image, ROSAT All-Sky Survey cutout, VLA FIRST 1.4 GHz, WISE Stern Donley AGN color cut, Eddington luminosity calculator, Shakura–Sunyaev disk temperature, fractional rms variability Vaughan 2003, Lomb–Scargle X-ray light curve, Little Red Dots JWST, L_X/L_UV Lusso Risaliti 2016, Reines & Volonteri 2015, hardness ratio AGN, Compton-thick Seyfert 2 NGC 1068, M87 jet, Cyg X-1 black hole binary, astroquery SkyView HEASARC, photutils aperture photometry, specutils Spectrum1D, FITS pitfalls.

---

## ⚡ Quick start

```bash
# 1. Clone and set up an environment (Python 3.11+)
git clone https://github.com/ysSemanticSystems/BlackHoleResearch.git
cd BlackHoleResearch
python3.11 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# 2. Download the curated multi-survey dataset (~50–200 MB from NASA SkyView)
python scripts/download_data.py

# 3. Launch the Streamlit app
streamlit run app.py
```

Browser opens at <http://localhost:8501>.

Want to add Chandra event files? See [Adding more data](#-adding-more-data) below for HEASARC, MAST, and IRSA recipes.

---

## 🔭 What you can do with this

### 🖼 Render any FITS image with WCS, properly

The right way: World Coordinate System overlay in J2000 sexagesimal, `origin='lower'` (FITS convention, not NumPy's default), and an **asinh stretch** so the dynamic range of an AGN core plus a galaxy halo plus an extended jet doesn't collapse to a black square. Supports `linear`, `sqrt`, `log`, `asinh`, and DS9-style `zscale`. Perceptually uniform colormaps only — no `jet`.

Backed by `astropy.visualization.ImageNormalize` and `WCSAxes`. The render path explicitly handles partial-or-missing WCS (a frequent issue in cutouts from older missions; see [PITFALLS #3](docs/PITFALLS.md)).

### 📡 X-ray event list-to-image binning

Chandra ACIS, XMM-Newton EPIC, NuSTAR FPM and NICER files are **photon tables**, not images. This app reads the `EVENTS` extension, picks the right `X`/`Y`/`ENERGY`(or `PI`) columns per mission, and bins them into a 2D image with an energy cut you control (soft <2 keV vs hard >2 keV by default — useful for spotting Compton-thick obscuration). The pitfalls of mistaking PI channels for eV (XMM/NuSTAR) are guarded against in code; see [PITFALLS #12](docs/PITFALLS.md).

### 📈 OGIP PHA spectrum loader with descriptive power-law fit

Loads OGIP-format `*.pha` / `*.pi` files, auto-detects `COUNTS` vs `RATE` columns, plots in log-log channel space with Poisson errors, and fits a power law `N(E) ∝ E^(−Γ)` over a channel range you select. **Honest disclaimer:** without RMF/ARF response files this is a *channel-space* slope, not a calibrated photon index — for publication-grade Γ, plug into Sherpa or XSPEC (planned Phase 2 hook). The classification table maps Γ ranges to AGN populations (Ricci+2017 BAT survey).

### 🌌 Multi-band SED builder: radio → mid-IR → optical → X-ray

The money plot for any accreting black hole. Takes points across nine decades of frequency (1.4 GHz to 50 keV), unit-checked via `astropy.units` (Jy ↔ erg/s/cm²/Hz ↔ νF_ν), color-coded by band, with an optional Elvis+1994 quasar template overlay. Resolves the radio jet synchrotron, the mid-IR torus bump, the optical/UV disk thermal continuum, and the X-ray corona power law on a single log-log axis. See `blackhole/sed.py`.

### ⏱ Light curve + Lomb–Scargle periodogram + Vaughan+2003 F_var

Bins arbitrary event lists into uniformly sampled light curves, computes a Lomb–Scargle periodogram via `astropy.timeseries.LombScargle` (the right tool for unevenly sampled, gap-filled data — see VanderPlas 2018), and computes the fractional rms variability amplitude with proper Vaughan+2003 Eq. B2 error bars. AGN typically show 5–30 % F_var on hour-day timescales; X-ray binaries can hit 50 % during bright hard states.

> ⚠️ **GTI handling is Phase 2 (M4).** Light curves over real Chandra/XMM/NuSTAR observations have *Good Time Interval* gaps that currently render as zero-count bins; full GTI masking is planned. See [PITFALLS #10](docs/PITFALLS.md).

### 🔬 Physics module with citations

Self-contained reference implementations:

| Module | Formulas | Primary source |
|---|---|---|
| [`accretion.py`](blackhole/physics/accretion.py) | Eddington luminosity, accretion rate Ṁ = L/(ηc²), Shakura–Sunyaev T(r) | Eddington 1926; Shakura & Sunyaev 1973 |
| [`spectral_xray.py`](blackhole/physics/spectral_xray.py) | Hardness ratio HR = (H−S)/(H+S) with Poisson error, photon-index classification | Ricci+2017 ApJS 233, 17 |
| [`infrared.py`](blackhole/physics/infrared.py) | Blackbody νB_ν(T), Stern+2012 single-color cut, Donley+2012 four-band wedge | Stern+2012, Donley+2012 |
| [`variability.py`](blackhole/physics/variability.py) | Excess variance, fractional rms F_var with Vaughan+2003 errors | Vaughan+2003 MNRAS 345, 1271 |

Every formula has its primary source in the docstring; enforcement of citation discipline is in [`.cursor/rules/01-documentation-standards.mdc`](.cursor/rules/01-documentation-standards.mdc).

### 📚 A documented catalog of common FITS failure modes

[`docs/PITFALLS.md`](docs/PITFALLS.md) is a 14-entry reference covering the recurring issues that come up when working with archival FITS files:

1. Event-list-as-image confusion
2. Linear stretch destroying dynamic range
3. WCS silently missing or partial
4. `origin='lower'` and flipped-image debugging
5. Cross-band unit confusion (Jy / mag / keV / Hz)
6. WISE Vega-vs-AB photometry
7. PHA spectra requiring RMF/ARF
8. Big-endian byte order in raw FITS arrays
9. Memmap invalidation after `with fits.open()` exits
10. GTI artifacts in light curves
11. Coordinate frame mismatches (ICRS / FK5 / Galactic)
12. PI-channel vs eV energy filters
13. Mission time systems (Chandra/XMM TT vs NuSTAR UTC)
14. Modified Julian Date conventions (the half-day offset)

Each entry has a symptom, a cause, and a concrete code defense with citations.

---

## 🎯 The targets and the physics

| Target | Why it's a milestone | Mass | Distance |
|---|---|---|---|
| **NGC 1068** | The archetypal **Seyfert 2 / Compton-thick AGN**. Antonucci & Miller 1985's polarized broad-line discovery here founded **AGN unification**. The mid-IR torus dominates everything < 30 µm. | M_BH ≈ 1.7 × 10⁷ M☉ (maser, Lodato+2003) | ~14 Mpc |
| **M87** | The first directly imaged **event horizon** (EHT Collaboration 2019). A low-luminosity AGN with an extended optical/radio jet observed since Curtis 1918. | M_BH = 6.5 × 10⁹ M☉ (EHT 2019) | ~16.8 Mpc |
| **Cyg X-1** | The first widely accepted **stellar-mass black hole** (Webster & Murdin 1972; Bolton 1972). Bright, variable, observed by every major X-ray mission since 1964. | M_BH ≈ 21 M☉ (Miller-Jones+2021) | ~2.2 kpc |

Each is downloaded as a stack of multi-band cutouts from NASA SkyView; see [`fits_data/MANIFEST.json`](fits_data/MANIFEST.json) for the SHA-stamped, provenance-tracked file list.

---

## 🌠 Phase 2 — JWST Little Red Dots research frontier

The post-2023 JWST surprise. **Little Red Dots (LRDs)** are compact, high-redshift sources with broad permitted lines implying accreting 10⁶–10⁸ M☉ black holes — but they show no detected X-ray emission (even in deep Chandra stacks; Akins+2025, Yue+2024) and no detected mid-IR torus. They appear to violate AGN unification in both bands at once.

Phase 2 (see [`PHASE2_PLAN.md`](PHASE2_PLAN.md)) adds:

- **LRD catalog** with peer-reviewed redshifts and M_BH (Greene+2024, Maiolino+2024, Harikane+2023, Akins+2025, Wang+2024).
- **JWST NIRSpec / NIRCam ingestion** via `astroquery.mast`.
- **Diagnostics** — exact reproductions of:
  - α_ox + L_X / L_UV vs **Lusso & Risaliti 2016 mean relation**;
  - **Balmer break detector** (Setton+2024, Wang+2024 definitions);
  - **Broad Hα → M_BH** via Reines & Volonteri 2015 calibration;
  - **Donley+2012 WISE wedge** with LRDs overlaid as the outliers they are.

The Phase 2 goal is a reproducible, citation-tracked LRD pipeline that emits calibrated values with uncertainties, and that anyone can run on a laptop.

---

## 🧬 Built on the astropy ecosystem

| Layer | Library | What we use it for |
|---|---|---|
| FITS I/O | [astropy.io.fits](https://docs.astropy.org/en/stable/io/fits/) | All file loading, header parsing, BinTable handling |
| Units | [astropy.units](https://docs.astropy.org/en/stable/units/) | End-to-end unit safety across SED bands |
| WCS | [astropy.wcs](https://docs.astropy.org/en/stable/wcs/) | World coordinate overlays on every image |
| Stretches | [astropy.visualization](https://docs.astropy.org/en/stable/visualization/) | Asinh, log, sqrt, ZScale + percentile clipping |
| Cosmology | [astropy.cosmology](https://docs.astropy.org/en/stable/cosmology/) (Phase 2) | Planck18 luminosity distances for LRDs |
| Periodograms | [astropy.timeseries](https://docs.astropy.org/en/stable/timeseries/) | Lomb–Scargle on unevenly sampled data |
| Archives | [astroquery](https://astroquery.readthedocs.io/) | SkyView, HEASARC, MAST, IRSA queries |
| Photometry | [photutils](https://photutils.readthedocs.io/) (Phase 2) | Aperture photometry + background annuli |
| Spectra | [specutils](https://specutils.readthedocs.io/) (Phase 2) | NIRSpec `_x1d` Spectrum1D handling |
| Sky regions | [regions](https://astropy-regions.readthedocs.io/) | DS9 region serialization |
| UI | [Streamlit](https://streamlit.io/) | Tabbed multi-target browser |

If you're an `astropy` user, you can read or extend this codebase fluently.

---

## 📂 Project structure

```text
BlackHoleResearch/
├── app.py                       # Streamlit UI — presentation layer only
├── PHASE2_PLAN.md               # Milestone roadmap with binary exit criteria
├── AGENTS.md                    # Contributor & AI-agent guide
├── CITATION.cff                 # Machine-readable citation metadata
├── codemeta.json                # Research software metadata (codemeta v2.0)
├── .zenodo.json                 # Zenodo upload metadata
├── README.md
├── LICENSE                      # MIT
├── requirements.txt
├── .cursor/rules/               # Persistent AI-agent guidance (10 focused rules)
├── .github/                     # (Phase 2 M0) CI workflows
├── fits_data/
│   ├── MANIFEST.json            # SHA-256, provenance, fetched-at per file
│   └── *.fits                   # Multi-survey cutouts (gitignored)
├── scripts/
│   └── download_data.py         # SkyView cross-mission cutout fetcher
├── blackhole/
│   ├── io.py                    # FITS load/inspect/event-list binning
│   ├── wcs_plot.py              # Image rendering w/ WCS + stretches
│   ├── spectra.py               # OGIP PHA loader + power-law fit
│   ├── sed.py                   # Multi-band SED assembly (units-safe)
│   ├── lightcurves.py           # Binning + Lomb-Scargle
│   └── physics/
│       ├── accretion.py         # Eddington, Ṁ, Shakura–Sunyaev
│       ├── spectral_xray.py     # Hardness ratio, photon-index classification
│       ├── infrared.py          # Blackbody, Stern+Donley WISE cuts
│       └── variability.py       # F_var, excess variance
└── docs/
    ├── PITFALLS.md              # 14-entry FITS failure-mode catalog
    └── adr/                     # Architecture Decision Records
```

---

## 🛰 Adding more data

### Chandra / XMM-Newton / NuSTAR (HEASARC)

```python
from astroquery.heasarc import Heasarc
import astropy.units as u

h = Heasarc()
obs = h.query_object("NGC 1068", catalog="chanmaster", radius=5 * u.arcmin)
print(obs[["obsid", "exposure", "ra", "dec"]])
h.download_data(obs[:1], host="heasarc", location="fits_data/")
```

`evt2.fits` event files drop directly into the **X-ray Events** tab.

### JWST / HST / TESS (MAST)

```python
from astroquery.mast import Observations

obs = Observations.query_object("CEERS-415", radius=0.05)
products = Observations.get_product_list(obs[:1])
Observations.download_products(products, download_dir="fits_data/")
```

NIRSpec `_x1d.fits` files will land in the Phase 2 NIRSpec viewer tab.

### Spitzer / WISE / 2MASS (IRSA)

```python
from astroquery.ipac.irsa import Irsa
result = Irsa.query_region("NGC 1068", catalog="allwise_p3as_psd",
                            radius="0d0m30s")
```

---

## 🤝 Contributing

Conventions are defined in [`AGENTS.md`](AGENTS.md) and enforced in CI.

Before opening a PR:

1. Read [`AGENTS.md`](AGENTS.md) for the project's prime directives and PR checklist.
2. Pick a milestone from [`PHASE2_PLAN.md`](PHASE2_PLAN.md) and tick its binary exit criteria in the PR description.
3. Run the local pre-merge checks (`ruff check`, `mypy --strict blackhole/`, `pytest -q`).

The persistent AI/contributor rules live in [`.cursor/rules/`](.cursor/rules/):

- `00-project-context.mdc` — repo identity (always loaded)
- `01-documentation-standards.mdc` — module/function docstring shape, citation requirements
- `02-scientific-honesty.mdc` — units, uncertainty, no invented values
- `03-python-style.mdc` — typing, dataclasses, error hierarchy
- `04-fits-handling.mdc` — FITS dos/donts mapped to PITFALLS
- `05-testing-discipline.mdc` — coverage targets, fixture conventions
- `06-ui-patterns.mdc` — Streamlit-as-presentation invariant
- `07-reproducibility.mdc` — provenance, manifest, caching, cosmology
- `08-adr.mdc` — when and how to write Architecture Decision Records
- `09-commit-and-pr.mdc` — commit message + PR template (always loaded)

---

## 📖 Citation

If you use this software in a publication, please cite via the **Cite this repository** widget on the GitHub sidebar (powered by [`CITATION.cff`](CITATION.cff)), or as:

```bibtex
@software{blackholeresearch_2026,
  author       = {ysSemanticSystems},
  title        = {{BlackHoleResearch: A Multiwavelength FITS Explorer for
                   Black Hole Astrophysics}},
  year         = {2026},
  url          = {https://github.com/ysSemanticSystems/BlackHoleResearch},
  version      = {0.1.0},
  license      = {MIT}
}
```

A Zenodo DOI is planned at the Phase 2 v0.2.0 release ([PHASE2_PLAN.md §M11](PHASE2_PLAN.md#m11--worked-notebooks--zenodo-release)). Please cite the **primary-source papers** for any physics formula you use (citations are inline in module docstrings and aggregated in [`PHASE2_PLAN.md §11`](PHASE2_PLAN.md#11-appendix--primary-source-references)).

---

## 🌐 See also

If this isn't the right tool for your task, one of these probably is:

- [**SAOImage DS9**](https://sites.google.com/cfa.harvard.edu/saoimageds9) — the canonical interactive FITS image viewer.
- [**CIAO + Sherpa**](https://cxc.cfa.harvard.edu/ciao/) — Chandra data reduction and forward-folded spectral fitting. The right tool for calibrated Γ measurements.
- [**XSPEC**](https://heasarc.gsfc.nasa.gov/xanadu/xspec/) — the X-ray spectral fitting standard, with response-file forward-folding.
- [**SAS**](https://www.cosmos.esa.int/web/xmm-newton/sas) — the XMM-Newton reduction package.
- [**JWST pipeline**](https://jwst-pipeline.readthedocs.io/) — official JWST data reduction (we consume its outputs, we don't replace it).
- [**Stingray**](https://stingray.science/) — Python X-ray timing analysis at depth (more advanced than our M4 light-curve module).
- [**aplpy**](https://aplpy.github.io/) — beautiful publication-quality WCS plots.
- [**Astropy tutorials**](https://learn.astropy.org/) — the canonical learning resource.
- [**NED**](https://ned.ipac.caltech.edu/) · [**SIMBAD**](https://simbad.u-strasbg.fr/simbad/) · [**ADS**](https://ui.adsabs.harvard.edu/) — the meta-archives every astronomer lives in.

This project's niche is **multiwavelength survey-data exploration with calibration-aware SEDs and LRD diagnostics**, with an emphasis on teaching the pitfalls.

---

## 📋 Standards adhered to

- **FITS Standard v4.0** — IAU FITS Working Group, 2018.
- **OGIP/92-007** — X-ray spectrum FITS conventions.
- **OGIP/93-003** — Type II PHA file structure.
- **WCS Paper II** (Calabretta & Greisen 2002) — World Coordinate System.
- **VOTable** — IVOA tabular data.
- **CFF 1.2.0** — Citation File Format.
- **CodeMeta 2.0** — Research Software Metadata.
- **SPDX MIT** — license identification.
- **Keep a Changelog** + **SemVer 2.0** — version & changelog conventions.

---

## 🪐 Acknowledgements

This software is built on the work of the **astropy** project (Astropy Collaboration 2013, 2018, 2022) and the **astroquery** project (Ginsburg et al. 2019, AJ 157, 98).

Data is provided by:

- **NASA SkyView** (HEASARC, Goddard Space Flight Center).
- **HEASARC** — High Energy Astrophysics Science Archive Research Center.
- **MAST** — Mikulski Archive for Space Telescopes (STScI).
- **IRSA** — NASA/IPAC Infrared Science Archive.
- **VizieR** — CDS, Strasbourg.

All NASA-archive data redistributed here is U.S. public domain.

---

<div align="center">

**Hosted on GitHub: <https://github.com/ysSemanticSystems/BlackHoleResearch>**

If this is useful to you, a ⭐ helps others find it.

</div>
