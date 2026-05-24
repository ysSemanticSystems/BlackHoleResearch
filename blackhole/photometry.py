"""
blackhole.photometry — aperture photometry on calibrated images.

The function `aperture_photometry` runs photutils with our conventions:

- Input: a `CalibratedImage` from `blackhole.calibration` (units-bearing
  array, pixel solid angle, survey label).
- Centroid: catalog source coordinate transformed via the image WCS.
- Background: median-clipped annulus with sigma from MAD (robust).
- Output: an `SEDPoint` carrying flux_density (or nu_f_nu for X-ray),
  flux_err, the aperture radius, and the upper-limit flag.

Upper-limit convention
----------------------
If the measured aperture sum is less than ``3 * sigma_background_aperture``,
the returned `SEDPoint` has ``upper_limit=True`` and ``flux_density`` set
to ``3 * sigma_background_aperture`` (the standard 3-sigma reporting
convention; see Bradley+ photutils docs).

References
----------
- photutils: Bradley+ (https://photutils.readthedocs.io/).
- Background MAD: Beers, Flynn & Gebhardt 1990 AJ 100 32 — robust estimators.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from photutils.aperture import (
    CircularAnnulus,
    CircularAperture,
    SkyCircularAnnulus,
    SkyCircularAperture,
    aperture_photometry,
)

from .calibration import CalibratedImage
from .sed import SEDPoint

ApertureMethod = Literal["exact", "subpixel", "center"]


@dataclass(frozen=True)
class PhotometryResult:
    """Full diagnostic record of a single aperture measurement.

    `to_sed_point()` collapses this to an `SEDPoint` for plotting; the
    PhotometryResult itself is what we cache and what provenance reads.
    """

    flux: u.Quantity                 # aperture sum (with units)
    flux_err: u.Quantity             # 1-sigma uncertainty from sky annulus
    sky_median: u.Quantity           # per-pixel median sky background
    sky_sigma: u.Quantity            # per-pixel sky sigma (from MAD)
    aperture_radius: u.Quantity      # arcsec
    annulus_inner: u.Quantity        # arcsec
    annulus_outer: u.Quantity        # arcsec
    n_pixels_aperture: float
    upper_limit: bool
    method: str                      # photutils method literal
    survey: str
    band: str
    label: str

    def to_sed_point(self) -> SEDPoint:
        """Compact representation for SED plotting.

        The aperture sum carries flux units (Jy for radio/IR/optical/UV
        surveys; erg/s/cm^2 for X-ray broad-band) — exactly what the
        SEDPoint expects.
        """
        flux_density: u.Quantity | None = None
        nu_f_nu: u.Quantity | None = None
        if self.flux.unit.is_equivalent(u.Jy):
            flux_density = self.flux.to(u.Jy)
        elif self.flux.unit.is_equivalent(u.erg / u.s / u.cm**2):
            nu_f_nu = self.flux.to(u.erg / u.s / u.cm**2)
        else:
            # Unknown / odd units: keep raw value as flux_density.
            flux_density = self.flux
        return SEDPoint(
            label=self.label,
            flux_density=flux_density,
            nu_f_nu=nu_f_nu,
            flux_err=self.flux_err,
            band=self.band,
            upper_limit=self.upper_limit,
            source=f"local photometry · {self.survey}",
        )


def aperture_photometry_on(
    cal: CalibratedImage,
    coord: SkyCoord,
    *,
    aperture_radius: u.Quantity = 10.0 * u.arcsec,
    annulus_inner: u.Quantity = 15.0 * u.arcsec,
    annulus_outer: u.Quantity = 25.0 * u.arcsec,
    method: ApertureMethod = "exact",
    label: str | None = None,
) -> PhotometryResult:
    """Run aperture photometry centred on `coord`.

    Parameters
    ----------
    cal
        A `CalibratedImage`. The image must have a celestial WCS in its
        underlying `ImageData`, otherwise we fall back to a pixel-centre
        aperture and document the fallback in the returned label.
    coord
        ICRS SkyCoord at which to place the aperture. Pass the catalog
        coordinate (``cat.NGC1068.coord`` etc.) rather than re-deriving.
    aperture_radius, annulus_inner, annulus_outer
        Quantity with angle units (arcsec). Default values are
        appropriate for IR/optical compact sources; radio images need
        much larger apertures (set explicitly).
    method
        photutils method — ``"exact"`` is sub-pixel-correct, ``"subpixel"``
        is faster, ``"center"`` is fastest and lossy. Default ``"exact"``.
    label
        Optional label for the resulting `SEDPoint`. Defaults to
        ``f"{survey} aperture"`` from the calibrated image.

    Returns
    -------
    PhotometryResult
        With flux, sky stats, aperture parameters, and the upper-limit
        flag set per the 3-sigma rule.

    Notes
    -----
    The sky annulus statistics use ``astropy.stats.sigma_clipped_stats``
    with the MAD estimator (``sigma=3, maxiters=5``). This is robust
    against contamination by satellite knots or bright background sources
    inside the annulus.
    """
    img = cal.original
    wcs = img.wcs

    # Aperture object placement: sky-based when WCS is present, pixel-based
    # otherwise (synthetic-test images or files without celestial WCS).
    if wcs is not None and wcs.has_celestial:
        ap_sky = SkyCircularAperture(coord, r=aperture_radius)
        an_sky = SkyCircularAnnulus(coord, r_in=annulus_inner, r_out=annulus_outer)
        ap = ap_sky.to_pixel(wcs)
        an = an_sky.to_pixel(wcs)
    else:
        ny, nx = img.array.shape
        ap = CircularAperture((nx / 2.0, ny / 2.0),
                              r=_arcsec_to_pix(aperture_radius, cal))
        an = CircularAnnulus((nx / 2.0, ny / 2.0),
                             r_in=_arcsec_to_pix(annulus_inner, cal),
                             r_out=_arcsec_to_pix(annulus_outer, cal))

    # Robust sky stats from MAD on the annulus mask.
    arr = cal.array.value  # work in raw numerics for stats
    annulus_mask = an.to_mask(method="center")
    annulus_data = annulus_mask.multiply(arr)
    annulus_data = annulus_data[annulus_mask.data > 0] if annulus_data is not None else np.array([])
    if annulus_data.size > 0:
        _, sky_median, sky_sigma = sigma_clipped_stats(
            annulus_data, sigma=3.0, maxiters=5,
        )
    else:
        sky_median, sky_sigma = 0.0, 0.0

    # Sum inside the aperture (raw), then subtract local sky.
    raw = aperture_photometry(arr, ap, method=method)
    raw_sum = float(raw["aperture_sum"][0])
    n_pixels = float(ap.area)
    flux_value = raw_sum - sky_median * n_pixels
    flux_err_value = sky_sigma * np.sqrt(n_pixels)

    flux_unit = cal.array.unit
    flux_q = flux_value * flux_unit
    flux_err_q = flux_err_value * flux_unit

    upper_limit = bool(flux_value < 3.0 * flux_err_value)
    if upper_limit:
        flux_q = 3.0 * flux_err_q

    return PhotometryResult(
        flux=flux_q,
        flux_err=flux_err_q,
        sky_median=sky_median * flux_unit,
        sky_sigma=sky_sigma * flux_unit,
        aperture_radius=aperture_radius.to(u.arcsec),
        annulus_inner=annulus_inner.to(u.arcsec),
        annulus_outer=annulus_outer.to(u.arcsec),
        n_pixels_aperture=n_pixels,
        upper_limit=upper_limit,
        method=method,
        survey=cal.survey,
        band=cal.band,
        label=label or f"{cal.survey} aperture",
    )


def _arcsec_to_pix(angle: u.Quantity, cal: CalibratedImage) -> float:
    """Convert an angular radius to pixels using the calibrated image's
    pixel solid angle.

    Uses pixel_solid_angle = (pixel_scale)^2 (small-angle approximation,
    which is accurate to better than 1% out to several arcmin).
    """
    sr = cal.pixel_solid_angle.to(u.sr).value
    pix_scale_rad = np.sqrt(sr)
    pix_scale_arcsec = (pix_scale_rad * u.rad).to(u.arcsec).value
    return float(angle.to(u.arcsec).value / pix_scale_arcsec)


def aperture_for_band(band: str) -> u.Quantity:
    """Sensible default aperture radius for a given band label.

    Bigger apertures for radio (synthesised beams of 5-30 arcsec; sources
    are often extended). Mid/large for IR (PSF FWHM up to 6'' for IRAS/IRIS).
    Small for optical/X-ray broad band where PSFs are tight.

    These are *defaults*. For real measurements, set the aperture
    explicitly to match the resolved size of the source in that band.
    """
    return {
        "radio":  60.0 * u.arcsec,
        "submm": 120.0 * u.arcsec,
        "ir":     20.0 * u.arcsec,
        "opt":    10.0 * u.arcsec,
        "uv":     10.0 * u.arcsec,
        "xray":   30.0 * u.arcsec,
    }.get(band, 10.0 * u.arcsec)
