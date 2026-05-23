# FITS Pitfalls Reference

A catalog of common failure modes when working with astronomical FITS data,
ordered roughly by how often they cause debugging time. Items 1–4 account
for the large majority of "why does my plot look wrong" issues.

---

## 1. X-ray data is not an image

Chandra, XMM-Newton, NuSTAR, and NICER all produce **event lists**: a
`BinTableHDU` where each row is one detected photon with columns `TIME`,
`X`, `Y`, `ENERGY` (or `PI`), and various status flags.

Calling `imshow` directly on this kind of HDU either errors or produces
meaningless output:

```python
hdul = fits.open("chandra_obs.fits")
plt.imshow(hdul[1].data)   # raises or returns garbage — that HDU is a Table, not an Image
```

To display an event list, bin the photon list into a 2D histogram first:

```python
events = bhio.load_events("chandra_obs.fits")
image, extent = bhio.bin_to_image(events, bins=512, energy_range_ev=(500, 7000))
plt.imshow(image, extent=extent, origin="lower")
```

The `blackhole.io` module handles this conversion. Most modern X-ray data
products arrive as event lists, so the distinction is worth internalizing.

---

## 2. Linear stretch destroys astronomical dynamic range

Astronomical images routinely span 4–7 orders of magnitude in pixel
intensity. A bright AGN core can be a million times brighter per pixel
than its host galaxy halo. Mapping that range linearly onto 256 display
levels collapses essentially all of it to either saturated white or pure
black.

The fix is a nonlinear stretch from `astropy.visualization`:

| Stretch | When to use |
|---|---|
| `LinearStretch` | When the data is already log-scaled (e.g. you took `log10` upstream) |
| `LogStretch` | Faint extended structure such as galaxy halos |
| `SqrtStretch` | General-purpose middle ground; common in HST press images |
| `AsinhStretch(a=0.1)` | Modern default; behaves linearly near zero, log-like at large values |
| `ZScaleInterval` + `LinearStretch` | The SAOImage DS9 default; reliable for quick inspection |

This module defaults to `asinh` because it handles both faint background
and bright cores without manual tuning.

---

## 3. WCS is silently missing or partial

Old plate-scan FITS files, small mission cutouts, and reduced data
products sometimes lack a complete World Coordinate System. The header
might have `CTYPE1='LINEAR'`, or no `CTYPE` keywords at all. Crucially,
`WCS(header)` returns a valid object regardless — it will just be a
meaningless identity mapping when the keywords are absent.

**Defense:** check that `wcs.has_celestial` is `True` and that
`wcs.naxis >= 2` before relying on RA/Dec. When in doubt, fall back to
pixel coordinates and label the axes honestly.

---

## 4. Image orientation: always pass `origin='lower'` to `imshow`

The FITS standard places pixel (1, 1) at the bottom-left of the image
on the sky. When astropy reads a FITS image into a NumPy array, the
array element `[0, 0]` corresponds to that bottom-left pixel. NumPy
arrays themselves have no inherent vertical orientation, but
`matplotlib.imshow` defaults to `origin='upper'`, which draws array
row 0 at the **top** of the figure.

The combination flips north–south unless you explicitly pass
`origin='lower'`:

```python
ax.imshow(data, origin='lower')   # correct: north is up
```

Forgetting this is the most common cause of "my source is on the wrong
side of the field" reports.

---

## 5. Unit confusion across wavebands

Different bands use different conventional flux units, and mixing them
in a single plot or calculation is a frequent source of orders-of-magnitude
errors.

| Band | Conventional flux unit |
|---|---|
| Radio / sub-mm | Jy (Jansky) = 10⁻²³ erg/s/cm²/Hz |
| IR / optical | Jy, mJy, or magnitudes (logarithmic, inverted scale) |
| UV | Jy or erg/s/cm²/Å |
| X-ray | erg/s/cm²/keV or count rate |
| Gamma-ray | photons/cm²/s/MeV |

**Defense:** use `astropy.units` consistently throughout. Convert to a
common form — typically νFν in erg/s/cm² — before plotting or comparing
across bands. The `sed.py` module enforces this through `astropy.units`
quantities at the point of ingestion.

---

## 6. WISE colors are reported in Vega, not AB

WISE catalogs report W1–W4 photometry in **Vega magnitudes**. The
Vega-to-AB offsets (Wright et al. 2010; WISE All-Sky Explanatory
Supplement, Section IV.4.h) are:

- W1 (3.4 μm): AB − Vega = 2.699
- W2 (4.6 μm): AB − Vega = 3.339
- W3 (12 μm):  AB − Vega = 5.174
- W4 (22 μm):  AB − Vega = 6.620

The Stern et al. (2012) and Donley et al. (2012) AGN color cuts are
defined in the Vega system. Applying them to AB-system photometry
without converting first will give incorrect classifications.

---

## 7. PHA spectra require response files for true energy axes

A PHA (Pulse Height Amplitude) spectrum file contains columns `CHANNEL`
and `COUNTS`. The channel index is a detector readout bin, not a
physical energy. Recovering true energies requires:

- **RMF** (Response Matrix File): the probabilistic mapping from channel
  to incident photon energy. One channel maps to a distribution of
  energies because of finite detector resolution.
- **ARF** (Auxiliary Response File): the effective area of the
  telescope + detector combination as a function of energy.

Production spectral fitting tools (XSPEC, Sherpa) forward-fold a
candidate model through the RMF and ARF and compare against observed
counts. Plotting raw channel against count rate and calling the result
"the spectrum" is a useful descriptive view, but it is not equivalent
to a calibrated spectrum.

**This app's power-law fits are intentionally descriptive.** A photon
index Γ of 1.9 fit in channel space is not the same as a Γ of 1.9 fit
through the proper response — the two will differ by mission and by
instrument. Use the values in this app for visualization and rough
classification, not for publication.

---

## 8. Big-endian byte order

FITS data is big-endian regardless of host architecture. Astropy handles
this transparently on access, so in normal usage there is nothing to do.
The issue surfaces only if you bypass the astropy API and pass the raw
data array through code that assumes native byte order (some older
C extensions, certain numpy operations on memory views).

**Defense:** stay within the astropy API for FITS data. If you must
operate on the raw array:

```python
arr = hdul[1].data.byteswap().newbyteorder()
```

---

## 9. Memory-map invalidation after the file is closed

`fits.open()` memory-maps array data by default. When you exit a `with`
block, the file closes and the memory map becomes invalid. Code that
holds a reference to `.data` outside the open block can segfault or
return corrupted values.

**Defense:** copy the array if it must outlive the open context:

```python
with fits.open(path) as hdul:
    data = hdul[0].data.copy()
# `data` is safe to use here
```

This module performs the copy in `load_image()` and `load_events()`.

---

## 10. Light-curve gaps are usually GTI artifacts, not source variability

Spacecraft observations are interrupted by orbit nights, Earth
occultations, instrument resets, South Atlantic Anomaly passages, and
similar non-source effects. The Good Time Interval (GTI) extension
records exactly when the instrument was accumulating valid data.
Binning event arrival times without applying the GTI produces
zero-count bins during the gaps, which look identical to real source
extinction.

**This app's `bin_events_to_lightcurve` does not currently apply GTI
filtering.** If a light curve appears to drop to zero counts for
thousands of seconds, the most likely explanation is a GTI gap rather
than source dimming. GTI integration is planned for Phase 3.

---

## 11. Coordinate frame mismatches

`SkyCoord("12h34m56s -01d23m45s", frame="icrs")` is not the same point
on the sky as the same string parsed as FK5 (different equinox handling),
and neither is the same as Galactic l, b. Cross-matching catalogs across
frames without explicit conversion produces silent positional errors at
the arcsecond level.

**Defense:** always specify the frame explicitly when constructing
`SkyCoord`. Convert between frames with `SkyCoord.transform_to(...)`
rather than re-parsing strings.

---

## 12. Energy filters on PI-channel data

Some missions store photon energies directly in eV (Chandra `energy`
column). Others store them as PI channel numbers (XMM-Newton, NuSTAR),
where the channel-to-energy conversion depends on the mission and
sometimes on event grade. Applying an `energy_range_ev=(500, 7000)`
filter to a PI-channel column silently filters on channel numbers
500–7000, which is not the same as the 0.5–7 keV energy band.

**Defense:** check the energy unit before filtering. The `bin_to_image`
function in this module applies eV filters only when
`events.energy_unit == "eV"`.

---

## 13. Time systems across missions

X-ray missions store photon arrival times as seconds since a
mission-specific reference epoch:

- **Chandra:** 1998-01-01 00:00:00 TT (Terrestrial Time)
- **XMM-Newton:** 1998-01-01 00:00:00 TT
- **NuSTAR:** 2010-01-01 00:00:00 UTC

To compare timestamps across missions, convert through `astropy.time.Time`
with the correct format and scale (`'tt'` versus `'utc'`). Naive
subtraction of raw mission-time values will give wrong absolute times,
including a 37-second offset between TT and UTC.

---

## 14. Modified Julian Date conventions

Many catalogs report observation epochs as Modified Julian Date:
MJD = JD − 2,400,000.5. The half-day offset matters: MJD starts at
midnight UTC, while JD starts at noon. Use astropy rather than
implementing the conversion by hand:

```python
from astropy.time import Time
t = Time(mjd_value, format='mjd', scale='utc')
```

---

## Adding to this list

This document is intended to grow as the project encounters new failure
modes. New entries should follow the existing structure: a short
description of the symptom, the underlying cause, a concrete defense or
code pattern, and a citation where relevant.
