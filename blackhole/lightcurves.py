"""
blackhole.lightcurves — time-series binning and periodogram analysis.

VARIABILITY IS THE MAIN HANDLE ON BH MASS AND SIZE
==================================================
Causality forces a black hole's variability timescale to be longer than
~R_g/c where R_g = GM/c². For a 10 M☉ stellar BH, R_g/c ≈ 50 µs — these
sources flicker on millisecond timescales (QPOs at tens of Hz).
For a 10⁸ M☉ supermassive BH, R_g/c ≈ 500 s — variability on hours to
years. Watching how a light curve flickers tells you about the mass, the
inner disk radius, and the corona geometry.

QUASI-PERIODIC OSCILLATIONS (QPOs)
==================================
Some X-ray binaries and a small number of AGN show narrow peaks in their
power spectrum — not strictly periodic but coherent enough that the QPO
frequency tracks the BH spin and inner disk radius via general relativistic
precession. The Lomb-Scargle periodogram (astropy.timeseries) handles
unevenly-sampled data, which X-ray data always is because of orbit gaps.

REFERENCES
==========
- van der Klis 2006, "Rapid X-ray Variability", in Compact Stellar X-ray
  Sources, eds. Lewin & van der Klis — the QPO bible.
- Vaughan et al. 2003, MNRAS 345, 1271 — fractional rms variability
  amplitude as the standard variability metric.
- VanderPlas 2018, ApJS 236, 16 — modern Lomb-Scargle for astronomy.
"""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import LombScargle
from matplotlib.figure import Figure

from .io import EventList


@dataclass
class LightCurve:
    times: np.ndarray   # seconds, relative to t_start
    rates: np.ndarray   # counts/sec per bin
    errors: np.ndarray  # 1-sigma
    bin_size_s: float
    energy_range_ev: tuple[float, float] | None
    source_label: str


def bin_events_to_lightcurve(
    events: EventList,
    bin_size_s: float = 100.0,
    energy_range_ev: tuple[float, float] | None = None,
) -> LightCurve:
    """Bin an event list into a uniform-time-step light curve.

    Each bin holds the photon count rate (counts/sec) and its Poisson
    error. Energy filtering is applied first if events.energy_unit is "eV".

    PITFALL: Event lists from real observations have *gaps* (spacecraft
    orbit nights, Earth occultations, instrument resets). Binning over
    gaps gives spurious zeros in your light curve. Production pipelines
    use Good Time Intervals (GTI) from a separate FITS extension to
    mask gaps; we ignore that for v1. Visible "rate drops to zero" in
    the resulting plot is usually a GTI artifact, not real source dimming.
    """
    t = events.times
    if energy_range_ev is not None and events.energies is not None and events.energy_unit == "eV":
        mask = (events.energies >= energy_range_ev[0]) & (events.energies <= energy_range_ev[1])
        t = t[mask]

    if t.size == 0:
        return LightCurve(
            times=np.array([]), rates=np.array([]), errors=np.array([]),
            bin_size_s=bin_size_s, energy_range_ev=energy_range_ev,
            source_label=f"{events.mission} (empty)")

    t_start = float(t.min())
    t_end = float(t.max())
    edges = np.arange(t_start, t_end + bin_size_s, bin_size_s)
    counts, _ = np.histogram(t, bins=edges)
    rates = counts / bin_size_s
    errors = np.sqrt(np.maximum(counts, 1.0)) / bin_size_s
    centers = 0.5 * (edges[:-1] + edges[1:]) - t_start

    return LightCurve(
        times=centers, rates=rates, errors=errors,
        bin_size_s=bin_size_s, energy_range_ev=energy_range_ev,
        source_label=events.mission,
    )


def render_lightcurve(
    lc: LightCurve,
    *,
    title: str | None = None,
    figsize: tuple[float, float] = (11, 5),
) -> Figure:
    """Render a light curve with error bars."""
    fig, ax = plt.subplots(figsize=figsize, facecolor="#0e1117")
    if lc.times.size == 0:
        ax.text(0.5, 0.5, "No events to plot", ha="center", va="center",
                color="white", transform=ax.transAxes, fontsize=14)
        return fig
    ax.errorbar(lc.times / 1000.0, lc.rates, yerr=lc.errors,
                fmt="o", ms=2, lw=0.5, color="#00f5d4", ecolor="#118ab2",
                alpha=0.7)
    ax.set_xlabel("Time since observation start (ks)", color="white", fontsize=12)
    ax.set_ylabel("Count rate (counts s$^{-1}$)", color="white", fontsize=12)
    band_str = ""
    if lc.energy_range_ev is not None:
        band_str = f"  {lc.energy_range_ev[0]/1000:.1f}-{lc.energy_range_ev[1]/1000:.1f} keV"
    ax.set_title(title or f"Light curve — {lc.source_label}{band_str}",
                 color="white", fontsize=13, pad=12)
    ax.set_facecolor("#0e1117")
    for spine in ax.spines.values():
        spine.set_color("white")
    ax.tick_params(colors="white")
    ax.grid(True, ls=":", alpha=0.3, color="white")
    fig.tight_layout()
    return fig


def lomb_scargle_periodogram(
    lc: LightCurve,
    min_freq_hz: float = 1e-5,
    max_freq_hz: float = 1e-1,
    n_freq: int = 5000,
) -> tuple[np.ndarray, np.ndarray]:
    """Lomb-Scargle periodogram of a light curve.

    Returns (frequencies_hz, power). Power is "normalized" Lomb-Scargle
    (peaks between 0 and 1, with 1 = perfect sinusoid).
    """
    if lc.times.size < 10:
        return np.array([]), np.array([])
    freqs = np.linspace(min_freq_hz, max_freq_hz, n_freq)
    ls = LombScargle(lc.times, lc.rates, dy=lc.errors)
    power = ls.power(freqs)
    return freqs, power


def render_periodogram(
    freqs: np.ndarray, power: np.ndarray,
    *,
    title: str | None = None,
    figsize: tuple[float, float] = (11, 5),
) -> Figure:
    """Render a Lomb-Scargle periodogram on log-log axes."""
    fig, ax = plt.subplots(figsize=figsize, facecolor="#0e1117")
    if freqs.size == 0:
        ax.text(0.5, 0.5, "Insufficient data for periodogram",
                ha="center", va="center", color="white",
                transform=ax.transAxes, fontsize=14)
        return fig
    ax.plot(freqs, power, color="#ff6b35", lw=1.0, alpha=0.9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Frequency (Hz)", color="white", fontsize=12)
    ax.set_ylabel("Lomb-Scargle power", color="white", fontsize=12)
    ax.set_title(title or "Lomb-Scargle periodogram",
                 color="white", fontsize=13, pad=12)
    ax.set_facecolor("#0e1117")
    for spine in ax.spines.values():
        spine.set_color("white")
    ax.tick_params(colors="white")
    ax.grid(True, ls=":", alpha=0.3, color="white", which="both")
    fig.tight_layout()
    return fig
