# BlackHoleResearch

> Multiwavelength FITS explorer for black hole astrophysics, designed for teaching and research workflows.

A Streamlit application that loads astronomical FITS files from multiple missions and wavebands, displays them with proper WCS overlays and intensity stretches, and assembles multi-band Spectral Energy Distributions from pixel-level data.

---

## What this is

Black holes do not emit light directly. They are detected through the radiation of accreting gas, which spans the entire electromagnetic spectrum:

| Band | Emitting region |
|---|---|
| Radio / sub-mm | Relativistic jets (synchrotron) |
| Infrared | Dusty torus surrounding the central engine |
| Optical / UV | Accretion disk (Shakura-Sunyaev thermal) |
| X-ray | Hot corona, disk reflection, jet base |
| Gamma-ray | Highest-energy jet particles |

No single band is sufficient on its own. This tool loads FITS files from multiple missions, displays each correctly, and assembles the **SED** that integrates them.

### Phase 1 (this build): classical regimes

- **NGC 1068** — archetypal Seyfert 2 / Compton-thick AGN
- **M87** — low-luminosity supermassive black hole with an extended optical/radio jet, and the first directly imaged event horizon (EHT 2019)
- **Cyg X-1** — first widely accepted stellar-mass black hole (Webster & Murdin 1972)

### Phase 2 (planned): Little Red Dots

The JWST research frontier — high-redshift sources with strong broad emission lines (implying accreting black holes) that show weak or absent X-ray and infrared emission, in tension with the standard AGN unification picture. References: Greene et al. 2024, Maiolino et al. 2024, Akins et al. 2025, Pacucci & Narayan 2025.

### Phase 3 (planned): calibrated spectral fitting

Sherpa or PyXspec integration with RMF/ARF response files for production-grade fits.

---

## Quick start

```bash
# 1. Set up an environment (Python 3.11 recommended)
python3.11 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# 2. Download the curated dataset (~50-200 MB from NASA SkyView)
python scripts/download_data.py

# 3. Launch the app
streamlit run app.py
```

Browser opens at `http://localhost:8501`.

---

## Project structure

```
BlackHoleResearch/
├── app.py                       # Streamlit UI
├── requirements.txt
├── README.md
├── .gitignore                   # FITS files are gitignored
├── fits_data/
│   ├── .gitkeep
│   └── MANIFEST.json            # records downloaded files
├── scripts/
│   └── download_data.py         # SkyView cross-mission cutouts
├── blackhole/
│   ├── __init__.py
│   ├── io.py                    # FITS loading + event-list binning
│   ├── wcs_plot.py              # image rendering w/ proper stretches
│   ├── spectra.py               # OGIP PHA spectrum loader + power-law fit
│   ├── sed.py                   # multi-band SED assembly
│   ├── lightcurves.py           # binning + Lomb-Scargle periodogram
│   └── physics/
│       ├── __init__.py
│       ├── accretion.py         # Eddington, M-dot, SS73 disk temperature
│       ├── spectral_xray.py     # hardness ratio, photon index classification
│       ├── infrared.py          # WISE color cuts (Stern, Donley), blackbody
│       └── variability.py       # F_var (Vaughan+2003)
└── docs/
    └── PITFALLS.md              # reference catalog of FITS handling issues
```

---

## The dataset

`scripts/download_data.py` pulls cutouts from [NASA SkyView](https://skyview.gsfc.nasa.gov/), HEASARC's cross-mission survey cutout service.

| Target | Surveys |
|---|---|
| NGC 1068 | DSS · 2MASS K · WISE W1/W3/W4 · ROSAT broad-band |
| M87 | DSS · 2MASS K · WISE W1 · VLA FIRST 1.4 GHz · ROSAT |
| Cyg X-1 | DSS · 2MASS K · WISE W1 · ROSAT |

Total download is approximately 50–200 MB depending on which surveys return data. The script is idempotent — repeated runs skip existing files — and `MANIFEST.json` records each download with a SHA-256 hash and the provenance URL.

For higher-resolution Chandra, XMM-Newton, and NuSTAR event files and OGIP spectra, see the **Adding more data** section below. These require direct HEASARC observation queries beyond the SkyView cutout interface.

---

## Adding more data

### From HEASARC (X-ray missions)

```python
from astroquery.heasarc import Heasarc
from astropy.coordinates import SkyCoord
import astropy.units as u

h = Heasarc()
# Find all Chandra observations of NGC 1068
obs = h.query_object("NGC 1068", catalog="chanmaster", radius=5*u.arcmin)
print(obs[["obsid", "exposure", "ra", "dec"]])

# Download a specific observation
h.download_data(obs[:1], host="heasarc", location="fits_data/")
```

The result includes `evt2` event files that load directly into the **X-ray Events** tab in the application.

### From MAST (HST, JWST, TESS)

```python
from astroquery.mast import Observations
obs = Observations.query_object("NGC 1068", radius=0.05)
products = Observations.get_product_list(obs[:1])
Observations.download_products(products, download_dir="fits_data/")
```

### From IRSA (Spitzer, WISE, 2MASS)

```python
from astroquery.ipac.irsa import Irsa
result = Irsa.query_region("NGC 1068", catalog="allwise_p3as_psd", radius="0d0m30s")
```

---

## The physics formulas in the code

All physics calculations have references back to original papers. See module docstrings for full citations. Key ones:

- **Eddington luminosity** (Eddington 1926): `L_Edd = 4πGMm_p·c / σ_T ≈ 1.26e38 (M/M☉) erg/s`
- **Shakura-Sunyaev disk** (1973): `T(r) ∝ (M·Ṁ)^¼ · r^(-¾)`
- **Hardness ratio**: `HR = (H-S)/(H+S)` — model-independent X-ray color
- **WISE AGN selection**: Stern+2012 single-color cut, Donley+2012 four-band wedge
- **Fractional rms variability** (Vaughan+2003): `F_var = √(σ²_xs)/⟨x⟩` with proper error bars

---

## Common pitfalls

See [`docs/PITFALLS.md`](docs/PITFALLS.md) for the full catalog. Top five:

1. **X-ray data is not an image.** Chandra/XMM/NuSTAR event files are tables of individual photons that must be binned into pixels before display.
2. **Linear stretch obscures the data.** Astronomical images span 4–7 orders of magnitude; a log, sqrt, or asinh stretch is required for meaningful display.
3. **WCS may be missing.** Old or cutout FITS files sometimes lack a complete World Coordinate System; verify `wcs.has_celestial` before assuming RA/Dec.
4. **PHA spectra require response files.** Channel space is not energy space. Power-law indices fit in channel space are descriptive, not calibrated.
5. **Light-curve gaps are typically GTI artifacts.** Real spacecraft observations have orbit gaps; the EVENTS extension alone does not record when the instrument was inactive.

---

## Why these targets

| Target | Significance |
|---|---|
| **NGC 1068** | The polarized broad-line discovery (Antonucci & Miller 1985) at this source established AGN unification. Its Compton-thick obscuration makes it a canonical case for the necessity of infrared diagnostics in identifying obscured AGN. |
| **M87** | First directly imaged event horizon (EHT Collaboration 2019). Hosts a ~1500-light-year optical/radio jet observed since Curtis 1918. A low-luminosity AGN, distinct from quasars but driven by the same fundamental accretion mechanism. |
| **Cyg X-1** | First widely accepted stellar-mass black hole (Webster & Murdin 1972; Bolton 1972). Bright, variable, and observed by every major X-ray mission since 1964. |

---

## License

MIT. NASA data is public domain.

## References

Primary references are inline in module docstrings. Key sources:

- FITS Standard v4.0 — IAU FITS Working Group, 2018
- Antonucci 1993, ARA&A 31, 473 — AGN unification
- Shakura & Sunyaev 1973, A&A 24, 337 — α-disk model
- Vaughan et al. 2003, MNRAS 345, 1271 — variability metrics
- Stern et al. 2012, ApJ 753, 30 — WISE AGN selection
- Donley et al. 2012, ApJ 748, 142 — refined WISE wedge
- Ricci et al. 2017, ApJS 233, 17 — modern AGN X-ray survey
- Pacucci & Narayan 2025 — super-Eddington LRD model
