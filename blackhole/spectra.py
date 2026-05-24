"""
blackhole.spectra — X-ray spectrum loaders and plots.

WHAT'S IN A "SPECTRUM" FROM AN X-RAY MISSION
============================================
An X-ray spectrum FITS file (PHA format, OGIP standard) is a BinTable HDU
with columns CHANNEL, COUNTS, STAT_ERR. The channels are detector readout
bins, NOT energies directly. To get energy you need the RMF (Response
Matrix File) which maps channel -> energy distribution. To get flux you
also need the ARF (Auxiliary Response File) which gives effective area
at each energy.

For v1 we plot in channel space and label honestly. Real spectral fitting
(power-law indices, blackbody fits, reflection features) requires XSPEC
or Sherpa with full response files — that's a Phase 3 deliverable. The
simple power-law fit here is purely descriptive.

OGIP STANDARD REFERENCES
========================
- George & Yusaf (1992), "The OGIP standard for the FITS format for
  spectra", OGIP/92-007:
  https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/ogip_92_007.html
- HEASARC OGIP memo summary:
  https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html

REFERENCES — POWER LAW FITS IN AGN
==================================
- Mushotzky, Done & Pounds (1993), "X-ray spectra and time variability
  of active galactic nuclei", ARA&A 31, 717 — the canonical review.
- For typical photon indices (Γ ≈ 1.7-2.1 for type-1 AGN, harder for
  obscured type-2), see e.g. Ricci et al. 2017, ApJS 233, 17 (BAT AGN
  Spectroscopic Survey).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib.figure import Figure
from scipy.optimize import curve_fit


@dataclass
class Spectrum:
    """A 1D X-ray spectrum (PHA-format)."""

    channels: np.ndarray
    counts: np.ndarray
    errors: np.ndarray | None  # 1-sigma per channel
    exposure_s: float | None
    mission: str
    instrument: str
    source_path: Path


def load_pha_spectrum(path: str | Path) -> Spectrum:
    """Load an OGIP-format PHA spectrum.

    Looks for the SPECTRUM extension first (OGIP standard), then any
    BinTable with a COUNTS or RATE column.

    PITFALL: Some spectra are stored as RATE (counts/sec) rather than
    raw COUNTS. We detect and load whichever is present. To convert
    rate -> counts you'd multiply by EXPOSURE; for plotting purposes
    either is fine, just label honestly.
    """
    path = Path(path)
    with fits.open(path) as hdul:
        spec_hdu = None
        for h in hdul:
            if h.name.upper() == "SPECTRUM" and isinstance(h, fits.BinTableHDU):
                spec_hdu = h
                break
        if spec_hdu is None:
            for h in hdul:
                if isinstance(h, fits.BinTableHDU):
                    names = [c.name.upper() for c in h.columns]
                    if "CHANNEL" in names and ("COUNTS" in names or "RATE" in names):
                        spec_hdu = h
                        break
        if spec_hdu is None:
            raise ValueError(f"No SPECTRUM extension in {path.name}")

        data = spec_hdu.data
        hdr = spec_hdu.header
        cols_upper = {c.name.upper(): c.name for c in spec_hdu.columns}

        channels = np.asarray(data[cols_upper["CHANNEL"]]).astype(float)
        if "COUNTS" in cols_upper:
            counts = np.asarray(data[cols_upper["COUNTS"]]).astype(float)
        else:
            counts = np.asarray(data[cols_upper["RATE"]]).astype(float)

        errors = None
        if "STAT_ERR" in cols_upper:
            errors = np.asarray(data[cols_upper["STAT_ERR"]]).astype(float)
        elif "COUNTS" in cols_upper:
            # Poisson approximation when no explicit error column.
            errors = np.sqrt(np.maximum(counts, 1.0))

        return Spectrum(
            channels=channels,
            counts=counts,
            errors=errors,
            exposure_s=hdr.get("EXPOSURE"),
            mission=(hdr.get("TELESCOP") or "UNKNOWN").upper(),
            instrument=(hdr.get("INSTRUME") or "UNKNOWN").upper(),
            source_path=path,
        )


def _power_law(channel: np.ndarray, norm: float, gamma: float) -> np.ndarray:
    """Simple power law N(E) ∝ E^(-Γ).

    We're fitting in channel space which is not strictly correct
    (channels are not linearly spaced in energy for most instruments),
    but for a *descriptive* photon-index estimate over a narrow energy
    band it's a reasonable first-pass. For publication-quality Γ you
    must use the RMF/ARF and a proper forward-folding fitter (XSPEC).
    """
    # +1 avoids divide-by-zero when channel == 0
    return norm * np.power(channel + 1.0, -gamma)


def fit_power_law(
    spectrum: Spectrum,
    channel_range: tuple[float, float] | None = None,
) -> tuple[float, float, np.ndarray]:
    """Fit a channel-space power law. Returns (alpha_channel, err, fit_curve).

    The returned scalar is the descriptive *channel-space* slope α_channel
    of N(channel) ∝ channel^(−α_channel), **not** a calibrated photon
    index Γ. Recovering Γ from PHA counts requires the RMF and ARF for
    the detector and a forward-folding fitter (XSPEC, Sherpa); see M5b
    in PHASE2_PLAN.md.

    Returns
    -------
    alpha_channel : float
        Best-fit channel-space slope.
    alpha_err : float
        1-sigma uncertainty (sqrt(pcov[1,1]) from curve_fit).
    fit_curve : np.ndarray
        Power-law evaluation over the full channel array.

    References
    ----------
    OGIP/92-007 — PHA file format; channel != energy in general.
    """
    ch = spectrum.channels
    co = spectrum.counts
    er = spectrum.errors if spectrum.errors is not None else np.sqrt(np.maximum(co, 1.0))

    if channel_range is not None:
        m = (ch >= channel_range[0]) & (ch <= channel_range[1]) & (co > 0)
    else:
        m = co > 0

    if m.sum() < 5:
        raise ValueError("Not enough non-zero channels to fit")

    p0 = (float(co[m].max()), 1.8)
    popt, pcov = curve_fit(_power_law, ch[m], co[m], p0=p0, sigma=er[m],
                           absolute_sigma=True, maxfev=5000)
    gamma = popt[1]
    gamma_err = float(np.sqrt(pcov[1, 1]))
    fit_curve = _power_law(ch, *popt)
    return gamma, gamma_err, fit_curve


def render_spectrum(
    spectrum: Spectrum,
    *,
    fit: tuple[float, float, np.ndarray] | None = None,
    title: str | None = None,
    figsize: tuple[float, float] = (10, 6),
) -> Figure:
    """Render a spectrum in log-log channel space with optional power-law fit."""
    fig, ax = plt.subplots(figsize=figsize, facecolor="#0e1117")
    m = spectrum.counts > 0

    if spectrum.errors is not None:
        ax.errorbar(spectrum.channels[m], spectrum.counts[m],
                    yerr=spectrum.errors[m], fmt="o", ms=2, lw=0.5,
                    color="#00d4ff", ecolor="#005f73", alpha=0.7,
                    label=f"{spectrum.instrument} data")
    else:
        ax.plot(spectrum.channels[m], spectrum.counts[m], "o", ms=2,
                color="#00d4ff", alpha=0.7, label=f"{spectrum.instrument} data")

    if fit is not None:
        alpha, alpha_err, fit_curve = fit
        ax.plot(spectrum.channels[m], fit_curve[m], "-", color="#ff6b35", lw=2,
                label=f"Power law: α_channel = {alpha:.2f} ± {alpha_err:.2f}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Channel (not energy)", color="white", fontsize=12)
    ax.set_ylabel("Counts" if spectrum.exposure_s else "Rate", color="white", fontsize=12)
    ax.set_title(title or f"{spectrum.mission} {spectrum.instrument} spectrum",
                 color="white", fontsize=13, pad=12)

    ax.legend(facecolor="#0e1117", edgecolor="white", labelcolor="white")
    ax.set_facecolor("#0e1117")
    for spine in ax.spines.values():
        spine.set_color("white")
    ax.tick_params(colors="white")
    ax.grid(True, ls=":", alpha=0.3, color="white")

    fig.tight_layout()
    return fig
