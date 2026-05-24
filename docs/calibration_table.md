# Calibration registry — supported surveys

`blackhole.calibration` flux-calibrates the surveys listed below. Anything
not in this table raises `UncalibratedDataError` rather than silently
returning a number.

Version: see `CALIBRATION_VERSION` in `blackhole/calibration.py`. Bump
when any constant changes.

| Survey key   | Input units      | Output units               | Required header keys      | Primary-source constants                              | Reference                                |
|--------------|------------------|----------------------------|---------------------------|-------------------------------------------------------|------------------------------------------|
| `2MASS-K`    | DN (data number) | Jy / pixel                 | `MAGZP`, CDELT or CD      | F_nu_0(K) = 666.7 Jy (Vega)                           | Cohen, Wheaton & Megeath 2003 AJ 126 1090 |
| `IRIS_12`    | MJy / sr         | Jy / pixel                 | `BUNIT='MJy/sr'`, CDELT   | — (no zero-point; per-pixel surface brightness)       | Miville-Deschênes & Lagache 2005 A&A 432 729 |
| `IRIS_25`    | MJy / sr         | Jy / pixel                 | `BUNIT='MJy/sr'`, CDELT   | —                                                     | Miville-Deschênes & Lagache 2005          |
| `IRIS_60`    | MJy / sr         | Jy / pixel (band='submm')  | `BUNIT='MJy/sr'`, CDELT   | —                                                     | Miville-Deschênes & Lagache 2005          |
| `IRIS_100`   | MJy / sr         | Jy / pixel (band='submm')  | `BUNIT='MJy/sr'`, CDELT   | —                                                     | Miville-Deschênes & Lagache 2005          |
| `AKARI`      | MJy / sr         | Jy / pixel                 | `BUNIT='MJy/sr'`, CDELT   | —                                                     | Doi+2015 PASJ 67 50                       |
| `VLA-FIRST`  | Jy / beam        | Jy / pixel                 | `BUNIT='Jy/beam'`, `BMAJ`, `BMIN`, CDELT | beam_area = π/(4 ln 2) · BMAJ · BMIN | Becker, White & Helfand 1995 ApJ 450 559  |
| `RASS`       | counts / pixel   | erg / s / cm² / pixel      | `EXPTIME`, CDELT          | ECF_broad = 1.08e-11 erg/s/cm² per (count/s)          | Snowden+1995 ApJ 454 643                  |
| `DSS`        | (raises)         | —                          | —                         | No absolute zero-point available                      | Lasker+1990 (DSS construction)            |

## Traps and caveats

- **IRIS pixel solid angle** is computed *per file* from CDELT1/CDELT2,
  not assumed. SkyView cutouts are usually 1.5 arcmin/pixel for IRIS but
  always verify against the header.
- **VLA FIRST beam** is approximately 5.4 arcsec FWHM but the actual
  BMAJ/BMIN are header-driven and slightly position-dependent. For
  unresolved sources, peak Jy/beam is the integrated flux. For extended
  emission, integrating Jy/pixel underestimates total flux unless the
  source is much smaller than the beam.
- **RASS ECF** here is the published broad-band value for a representative
  T ≈ 1 keV thermal spectrum with N_H ≈ 3×10²⁰ cm⁻². The true ECF
  depends on the source spectrum and absorbing column. For Galactic XRBs
  with high N_H (Cyg X-1: N_H ≈ 6×10²¹ cm⁻²) the broad-band ECF can be
  off by tens of percent; use Snowden+1995 Fig. 8 if you need a better
  value, or fit XSPEC/Sherpa with a real RMF.
- **2MASS MAGZP** is per-image and lives in the FITS header for SkyView
  cutouts. If you obtained the image through another route, verify
  MAGZP is present before calling the calibrator.
- **DSS** raises because absolute photometric calibration of plate scans
  requires per-plate response correction against catalog stars, which is
  a research task not a one-liner. Use DSS for visual context only.

## Calibration-version policy

`CALIBRATION_VERSION` is the contract that downstream code (photometry,
provenance, cached results) checks. Bump rules:

- **Patch** (1.0.0 → 1.0.1): bug fixes that do not change values for
  a previously calibrated image (e.g. a docstring fix, a clearer error
  message).
- **Minor** (1.0.0 → 1.1.0): new survey added, no existing values changed.
- **Major** (1.0.0 → 2.0.0): a constant changed (F_nu_0, ECF) or a
  formula was corrected. All caches keyed on version are invalidated.
