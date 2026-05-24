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

from dataclasses import dataclass, field

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
    effective_exposure: np.ndarray = field(default_factory=lambda: np.array([]))
    gti_applied: bool = False
    total_exposure_s: float = 0.0  # sum of effective_exposure across all bins
    gti_total_s: float = 0.0       # total GTI duration over the binned span


def _gti_overlap_per_bin(
    bin_edges: np.ndarray, gti: np.ndarray,
) -> np.ndarray:
    """Total exposure inside each bin, intersected with the GTIs.

    Parameters
    ----------
    bin_edges
        Strictly increasing array of length N+1 defining N bins.
    gti
        (M, 2) array of (start, stop) intervals, sorted and non-overlapping.

    Returns
    -------
    np.ndarray
        Length-N array of seconds of effective exposure per bin.
    """
    n_bins = bin_edges.size - 1
    out = np.zeros(n_bins, dtype=float)
    if gti.size == 0:
        return out
    gti_starts = gti[:, 0]
    gti_stops = gti[:, 1]
    for i in range(n_bins):
        bs, be = bin_edges[i], bin_edges[i + 1]
        if be <= gti_starts[0] or bs >= gti_stops[-1]:
            continue
        lo = np.searchsorted(gti_stops, bs, side="right")
        hi = np.searchsorted(gti_starts, be, side="left")
        if lo >= hi:
            continue
        starts = np.maximum(gti_starts[lo:hi], bs)
        stops = np.minimum(gti_stops[lo:hi], be)
        seg = np.clip(stops - starts, 0.0, None)
        out[i] = float(seg.sum())
    return out


def bin_events_to_lightcurve(
    events: EventList,
    bin_size_s: float = 100.0,
    energy_range_ev: tuple[float, float] | None = None,
    *,
    apply_gti: bool = True,
    min_exposure_fraction: float = 0.5,
) -> LightCurve:
    """Bin an event list into a uniform-time-step light curve.

    Each bin holds the photon count rate (counts/sec) and its Poisson
    error. Energy filtering is applied first if events.energy_unit is "eV".

    Good Time Interval handling
    ---------------------------
    If ``apply_gti=True`` and ``events.gti`` is populated, the per-bin
    effective exposure is computed as the intersection of the bin with
    the GTIs. Bins whose effective exposure / bin_size < ``min_exposure_fraction``
    are *dropped from the output* (not rendered as zero). The count rate
    and Poisson error of surviving bins are normalized by the effective
    exposure, not by the nominal bin width.

    Without GTIs ``apply_gti=False`` or ``events.gti is None``) the
    effective exposure is set to ``bin_size_s`` for every bin and the
    historical behavior is preserved.

    References
    ----------
    Vaughan+2003 MNRAS 345 1271 — exposure-weighted error treatment.
    Chandra Proposers' Observatory Guide §4 — GTI semantics.
    """
    t = events.times
    if energy_range_ev is not None and events.energies is not None and events.energy_unit == "eV":
        mask = (events.energies >= energy_range_ev[0]) & (events.energies <= energy_range_ev[1])
        t = t[mask]

    if t.size == 0:
        return LightCurve(
            times=np.array([]), rates=np.array([]), errors=np.array([]),
            bin_size_s=bin_size_s, energy_range_ev=energy_range_ev,
            source_label=f"{events.mission} (empty)",
            effective_exposure=np.array([]),
            gti_applied=False,
            total_exposure_s=0.0,
            gti_total_s=0.0,
        )

    gti = events.gti if apply_gti else None
    use_gti = gti is not None and gti.size > 0
    if use_gti and gti is not None:
        # Span the full GTI window so we don't miss exposure at the edges
        # (the first/last event may not land near the GTI boundaries).
        t_start = float(min(t.min(), float(gti[0, 0])))
        t_end = float(max(t.max(), float(gti[-1, 1])))
    else:
        t_start = float(t.min())
        t_end = float(t.max())
    edges = np.arange(t_start, t_end + bin_size_s, bin_size_s)
    counts, _ = np.histogram(t, bins=edges)

    if use_gti and gti is not None:
        gti_total = float((gti[:, 1] - gti[:, 0]).sum())
        effective = _gti_overlap_per_bin(edges, gti)
        keep = effective / bin_size_s >= min_exposure_fraction
        if not np.any(keep):
            return LightCurve(
                times=np.array([]), rates=np.array([]), errors=np.array([]),
                bin_size_s=bin_size_s, energy_range_ev=energy_range_ev,
                source_label=f"{events.mission} (no GTI-overlapping bins)",
                effective_exposure=np.array([]),
                gti_applied=True,
                total_exposure_s=0.0,
                gti_total_s=gti_total,
            )
        counts = counts[keep]
        effective = effective[keep]
        centers = 0.5 * (edges[:-1] + edges[1:]) - t_start
        centers = centers[keep]
        rates = counts / effective
        errors = np.sqrt(np.maximum(counts, 1.0)) / effective
        return LightCurve(
            times=centers, rates=rates, errors=errors,
            bin_size_s=bin_size_s, energy_range_ev=energy_range_ev,
            source_label=events.mission,
            effective_exposure=effective,
            gti_applied=True,
            total_exposure_s=float(effective.sum()),
            gti_total_s=gti_total,
        )

    rates = counts / bin_size_s
    errors = np.sqrt(np.maximum(counts, 1.0)) / bin_size_s
    centers = 0.5 * (edges[:-1] + edges[1:]) - t_start
    effective = np.full(counts.shape, bin_size_s, dtype=float)
    return LightCurve(
        times=centers, rates=rates, errors=errors,
        bin_size_s=bin_size_s, energy_range_ev=energy_range_ev,
        source_label=events.mission,
        effective_exposure=effective,
        gti_applied=False,
        total_exposure_s=float(effective.sum()),
        gti_total_s=0.0,
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
    gti_str = ""
    if lc.gti_applied and lc.gti_total_s > 0:
        frac = lc.total_exposure_s / lc.gti_total_s
        gti_str = f"  ·  GTI exposure {lc.total_exposure_s:.0f}s ({frac*100:.0f}% of GTI window)"
    ax.set_title(title or f"Light curve — {lc.source_label}{band_str}{gti_str}",
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
