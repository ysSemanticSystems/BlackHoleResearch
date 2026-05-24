"""
blackhole.sed — Spectral Energy Distribution assembly across wavebands.

WHAT AN SED IS AND WHY IT'S THE MONEY PLOT
==========================================
A Spectral Energy Distribution shows flux (or νF_ν) as a function of
frequency / wavelength / energy across the entire electromagnetic spectrum
on log-log axes. For black holes this is the single most diagnostic plot
because different physical components peak in different bands:

  Radio (GHz):           Synchrotron from jets
  Far-IR (100 µm):       Cool dust in host galaxy
  Mid-IR (10 µm):        Hot dust torus (the AGN signature)
  Near-IR (1-3 µm):      Hottest dust + stellar light
  Optical/UV:            Accretion disk (Shakura-Sunyaev, T ~ 10^4-10^5 K)
  Soft X-ray (< 2 keV):  Disk inner edge / soft excess
  Hard X-ray (> 2 keV):  Corona (inverse Compton)
  Gamma-ray:             Jet inverse Compton, pair processes

If the SED has a "big blue bump" in UV plus a hard X-ray power law plus a
mid-IR torus bump, it's a classical AGN. If the IR and X-ray are both weak
relative to the optical, it's something weirder — like a Little Red Dot.

UNITS — THE NEW USER FOOTGUN
============================
Different bands report flux in different conventional units:

  Radio/sub-mm: Jansky (Jy = 10^-23 erg/s/cm²/Hz)
  IR/optical:  Jy, mJy, or mag (logarithmic, inverted!)
  UV/X-ray:    erg/s/cm²/Hz or erg/s/cm²/keV
  Gamma-ray:   photons/cm²/s/MeV

The astropy.units system enforces conversion. Convert everything to a
common νF_ν in erg/s/cm² before plotting. That's what `to_nu_fnu` does.

REFERENCES
==========
- Elvis et al. 1994, "Atlas of Quasar Energy Distributions", ApJS 95, 1.
  The canonical broad-band quasar SED template.
- Shang et al. 2011, ApJS 196, 2 — modernized SED templates.
- Stern et al. 2012, ApJ 753, 30 — WISE mid-IR AGN color selection.
- Donley et al. 2012, ApJ 748, 142 — refined WISE wedge.
- Lusso & Risaliti 2016, ApJ 819, 154 — the L_X / L_UV relation that
  Little Red Dots violate.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import astropy.units as u
import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as const
from matplotlib.figure import Figure


@dataclass
class SEDPoint:
    """One photometric measurement to drop on an SED.

    Provide ONE of (wavelength, frequency, energy). The others get filled.
    Provide flux as flux_density (per Hz) or already-converted nu_f_nu.
    """

    label: str                       # "WISE W1", "Chandra 2-10 keV", etc.
    wavelength: u.Quantity | None = None   # e.g. 3.4 * u.micron
    frequency: u.Quantity | None = None
    energy: u.Quantity | None = None       # for X-rays: 5 * u.keV
    flux_density: u.Quantity | None = None  # Jy or erg/s/cm²/Hz
    nu_f_nu: u.Quantity | None = None       # erg/s/cm²  (already converted)
    flux_err: u.Quantity | None = None
    band: str = "unknown"            # "radio", "ir", "opt", "uv", "xray", "gamma"
    upper_limit: bool = False        # True = plot as down-arrow
    source: str = ""                 # provenance: archive, paper, etc.


@dataclass
class SED:
    """A collection of SEDPoints for one astronomical source."""

    target_name: str
    points: list[SEDPoint] = field(default_factory=list)
    redshift: float | None = None    # for cosmological sources

    def add(self, point: SEDPoint) -> SED:
        self.points.append(point)
        return self


def to_frequency(p: SEDPoint) -> u.Quantity:
    """Resolve a point's frequency from whichever of (wavelength, frequency,
    energy) was supplied."""
    if p.frequency is not None:
        return p.frequency.to(u.Hz)
    if p.wavelength is not None:
        return (const.c / p.wavelength).to(u.Hz)
    if p.energy is not None:
        # E = h ν
        return (p.energy / const.h).to(u.Hz)
    raise ValueError(f"SEDPoint {p.label}: need wavelength, frequency, or energy")


def to_nu_fnu(p: SEDPoint) -> u.Quantity:
    """Convert to νF_ν in erg/s/cm². If already given, return it."""
    if p.nu_f_nu is not None:
        return p.nu_f_nu.to(u.erg / u.s / u.cm**2)
    if p.flux_density is None:
        raise ValueError(f"SEDPoint {p.label}: need flux_density or nu_f_nu")
    nu = to_frequency(p)
    fnu = p.flux_density.to(u.erg / u.s / u.cm**2 / u.Hz)
    return (nu * fnu).to(u.erg / u.s / u.cm**2)


# Band -> color (matched to a dark-mode aesthetic)
BAND_COLORS = {
    "radio": "#ff006e",
    "submm": "#8338ec",
    "ir":    "#fb5607",
    "opt":   "#ffbe0b",
    "uv":    "#3a86ff",
    "xray":  "#00f5d4",
    "gamma": "#f72585",
    "unknown": "#aaaaaa",
}


def render_sed(
    sed: SED,
    *,
    title: str | None = None,
    figsize: tuple[float, float] = (11, 7),
    show_legend: bool = True,
    overplot_quasar_template: bool = False,
) -> Figure:
    """Render the SED with band-coded markers on log-log axes.

    If `overplot_quasar_template=True`, draws an Elvis et al. 1994 mean
    radio-loud quasar SED for visual comparison. (Hand-tabulated template
    values; for publication overlay a proper template from VizieR.)
    """
    fig, ax = plt.subplots(figsize=figsize, facecolor="#0e1117")

    if not sed.points:
        ax.text(0.5, 0.5, "No SED points loaded",
                ha="center", va="center", color="white",
                transform=ax.transAxes, fontsize=14)
        return fig

    # Group by band so the legend stays tidy.
    by_band: dict[str, list[SEDPoint]] = {}
    for p in sed.points:
        by_band.setdefault(p.band, []).append(p)

    for band, pts in by_band.items():
        nus = np.array([to_frequency(p).value for p in pts])
        nufnu = np.array([to_nu_fnu(p).value for p in pts])
        upper = np.array([p.upper_limit for p in pts])
        color = BAND_COLORS.get(band, "#aaaaaa")

        # Detected points
        if (~upper).any():
            ax.scatter(nus[~upper], nufnu[~upper], s=80,
                       color=color, edgecolor="white", linewidth=0.5,
                       label=band.upper(), zorder=3)
        # Upper limits as down arrows
        if upper.any():
            ax.scatter(nus[upper], nufnu[upper], s=120,
                       color=color, marker="v", edgecolor="white",
                       linewidth=0.5, label=f"{band.upper()} (upper limit)",
                       zorder=3)

    if overplot_quasar_template:
        _overplot_elvis_template(ax)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Frequency ν (Hz)", color="white", fontsize=12)
    ax.set_ylabel(r"$\nu F_\nu$ (erg s$^{-1}$ cm$^{-2}$)", color="white", fontsize=12)

    full_title = title or f"{sed.target_name} — Spectral Energy Distribution"
    if sed.redshift is not None:
        full_title += f"  (z = {sed.redshift:.3f})"
    ax.set_title(full_title, color="white", fontsize=14, pad=14)

    # Secondary axis: wavelength. Matplotlib calls these with array-like
    # and we return numpy arrays. The stub's Callable signature is too
    # narrow for our np.asarray-coerced inputs, so type checks are
    # suppressed locally.
    def freq_to_wave_um(nu_hz: np.ndarray) -> np.ndarray:
        with np.errstate(divide="ignore"):
            return np.asarray((const.c.value / np.asarray(nu_hz)) * 1e6)

    def wave_um_to_freq(lam_um: np.ndarray) -> np.ndarray:
        with np.errstate(divide="ignore"):
            return np.asarray(const.c.value / (np.asarray(lam_um) * 1e-6))

    secax = ax.secondary_xaxis(
        "top",
        functions=(freq_to_wave_um, wave_um_to_freq),  # type: ignore[arg-type]
    )
    secax.set_xlabel("Wavelength (µm)", color="white", fontsize=11)
    secax.tick_params(colors="white")

    if show_legend:
        ax.legend(facecolor="#0e1117", edgecolor="white",
                  labelcolor="white", loc="best")

    ax.set_facecolor("#0e1117")
    for spine in ax.spines.values():
        spine.set_color("white")
    ax.tick_params(colors="white")
    ax.grid(True, ls=":", alpha=0.3, color="white")
    fig.tight_layout()
    return fig


def _overplot_elvis_template(ax: matplotlib.axes.Axes) -> None:
    """Schematic Elvis et al. 1994 mean RL-quasar νF_ν template (normalized).

    Hand-tabulated landmarks (log νFν vs log ν, arbitrary normalization).
    For real overlay work, pull the digitized template from VizieR catalog
    J/ApJS/95/1.
    """
    log_nu = np.array([9, 11, 13, 14, 15, 16, 17, 18, 19])     # log10(Hz)
    log_nufnu = np.array([-2.5, -1.8, -1.2, -1.0, -0.6, -0.8, -1.2, -1.5, -2.0])
    # Renormalize: shift to a useful position relative to actual data range
    nu = 10**log_nu
    nufnu = 10**log_nufnu * 1e-10  # arbitrary visual scale
    ax.plot(nu, nufnu, "--", color="white", alpha=0.4, lw=1.5,
            label="Elvis+1994 quasar (schematic)")
