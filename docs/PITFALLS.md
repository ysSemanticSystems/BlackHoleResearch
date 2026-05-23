# FITS Pitfalls — the file you'll re-read at 2am

Every trap that bites people working with astronomical FITS data, in roughly the order they bite. The first three are responsible for ~80% of "why is my plot broken" debugging.

---

## 1. X-ray "data" is not an image

Chandra, XMM-Newton, NuSTAR, NICER all produce **event lists**: a `BinTableHDU` where each row is one detected photon with columns `TIME`, `X`, `Y`, `ENERGY` (or `PI`), and various flags.

**Wrong:**
```python
hdul = fits.open("chandra_obs.fits")
plt.imshow(hdul[1].data)   # raises or shows garbage — that's a Table, not an Image
```

**Right:**
```python
events = bhio.load_events("chandra_obs.fits")
image, extent = bhio.bin_to_image(events, bins=512, energy_range_ev=(500, 7000))
plt.imshow(image, extent=extent, origin="lower")
```

This module handles it for you, but you'll meet raw FITS often enough that knowing the convention matters.

---

## 2. Linear stretch turns astronomy into a black square

Astronomical images have 4-7 orders of magnitude dynamic range. A bright AGN core can be a million times brighter per pixel than its host galaxy halo. Linearly mapping that to 256 display levels collapses everything to either pure white or pure black.

**Fix:** use a nonlinear stretch via `astropy.visualization`:

| Stretch | When to use |
|---|---|
| `LinearStretch` | When the data is already log (e.g. you took log10 yourself) |
| `LogStretch` | Faint structure, galaxy halos |
| `SqrtStretch` | Middle ground, Hubble press-release default |
| `AsinhStretch(a=0.1)` | Modern default; well-behaved near zero, log-like at large values |
| `ZScaleInterval + LinearStretch` | The DS9 default; usually works |

This module defaults to `asinh` because it's the most forgiving across data types.

---

## 3. WCS is silently missing or partial

Old plate-scan FITS files, small mission cutouts, and reduced data products sometimes lack World Coordinate System keywords. The header might have `CTYPE1='LINEAR'` or no `CTYPE` at all. `WCS(header)` returns a valid object regardless — but it'll be a meaningless identity mapping.

**Defense:** check `wcs.naxis >= 2` and verify the projection looks sensible. If unsure, fall back to pixel coordinates and label the plot honestly.

---

## 4. Origin convention

FITS pixel `(1,1)` is the **bottom-left** corner (Fortran convention). numpy arrays are indexed from `(0,0)` at the **top-left**. matplotlib `imshow` defaults to numpy convention.

**Always pass `origin='lower'`** to `imshow` when displaying FITS images. Forgetting flips north-south.

---

## 5. Unit confusion across bands

| Band | Conventional flux unit |
|---|---|
| Radio / sub-mm | Jy (Jansky) = 10⁻²³ erg/s/cm²/Hz |
| IR / optical | Jy, mJy, or magnitudes (logarithmic, inverted) |
| UV | Jy or erg/s/cm²/Å |
| X-ray | erg/s/cm²/keV or count rate |
| Gamma-ray | photons/cm²/s/MeV |

**Defense:** use `astropy.units` consistently. Convert to a common form (`νF_ν` in erg/s/cm²) before plotting across bands. The `sed.py` module enforces this.

### 5a. WISE colors are Vega, not AB

WISE catalogs report W1-W4 photometry in **Vega magnitudes**. AB conversions (Wright et al. 2010):
- W1: AB - Vega = 2.699
- W2: AB - Vega = 3.339
- W3: AB - Vega = 5.174
- W4: AB - Vega = 6.620

The Stern and Donley AGN cuts are defined in Vega — don't apply them to AB-magnitude photometry without converting first.

---

## 6. PHA spectra need response files

A PHA (Pulse Height Amplitude) spectrum file has columns `CHANNEL` and `COUNTS`. The channel is a detector readout bin, not an energy. To get true energies you need:

- **RMF** (Response Matrix File): channel → energy distribution (one channel maps to multiple energies with different probabilities due to detector response)
- **ARF** (Auxiliary Response File): effective area vs energy (telescope+detector throughput)

Real spectral fitting (XSPEC, Sherpa) forward-folds a model through RMF and ARF. Plotting raw channels vs counts and calling that "the spectrum" is descriptive only.

**The power-law fits in this app are deliberately descriptive.** A Γ of 1.9 fit from channels is not the same as a Γ of 1.9 fit through an RMF — they differ by mission and instrument. Use those values for visualization context only.

---

## 7. Big-endian byte order

FITS data is big-endian regardless of host machine. astropy handles this transparently on access, but if you save a numpy array to disk in mixed-architecture code paths you can get silent corruption.

**Defense:** stay in the astropy API. If you must do raw numpy on FITS arrays:
```python
arr = hdul[1].data.byteswap().newbyteorder()
```

---

## 8. Memory map invalidation

`fits.open()` memory-maps data by default. When you exit the `with` block, the file closes and the memmap is invalid. Code that accesses `.data` outside the block can segfault or return garbage.

**Defense:** always `.copy()` the array if you'll use it outside the open block:
```python
with fits.open(path) as hdul:
    data = hdul[0].data.copy()
# `data` is safe here
```

This module does this for you in `load_image()` and `load_events()`.

---

## 9. Light-curve gaps are GTI artifacts, not source variability

Spacecraft observations have orbit nights, Earth occultations, instrument resets, SAA passages. The Good Time Interval (GTI) extension tells you when the instrument was actually accumulating data. Binning event arrival times without applying the GTI mask gives you spurious zero-rate bins in the gaps.

**This app's `bin_events_to_lightcurve` does not apply GTI** — for v1 simplicity. If you see "the source dropped to zero counts for 4000 seconds," that's almost certainly a GTI gap, not real source dimming. Phase 3 will integrate the GTI.

---

## 10. Coordinate frame mismatches

`SkyCoord("12h34m56s -01d23m45s", frame="icrs")` is not the same coordinate as the same string in `frame="fk5"` (drift since the equinox), and definitely not the same as Galactic l/b. Cross-matching catalogs without checking frames produces silent ~1 arcsec errors.

**Defense:** be explicit about frames; convert via `SkyCoord.transform_to(...)`.

---

## 11. Energy filters on PI-channel data

Some missions store energies in eV (Chandra) and others in PI channels (XMM, NuSTAR — where the conversion to eV is mission- and grade-dependent). Applying an `energy_range_ev=(500, 7000)` filter to PI-channel data silently filters on channel numbers 500-7000, which is *not* the same energy range.

**Defense:** check `events.energy_unit` before filtering. The `bin_to_image` function only applies eV filters when `energy_unit == "eV"`.

---

## 12. Naive vs astropy time

X-ray missions store times as seconds since a mission-specific epoch (Chandra: 1998-01-01 TT; XMM: 1998-01-01 TT; NuSTAR: 2010-01-01 UTC). To compare across missions, convert via `astropy.time.Time` with the correct format and scale. Naive subtraction gives wrong absolute times.

---

## 13. The MJD obsession

Many old catalogs report observation epochs as Modified Julian Date. MJD = JD - 2400000.5. The 2400000.5 offset is a convention — half-day offsets exist (MJD starts at midnight UTC, JD at noon). Use `Time(mjd_value, format='mjd', scale='utc')`, don't roll your own.

---

If you hit something not in this list, add it. The list is the documentation.
