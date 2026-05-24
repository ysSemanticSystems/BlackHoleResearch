"""
blackhole.wcs_plot — image rendering with proper astronomical conventions.

WHY THIS MODULE EXISTS
======================
Astronomical images have dynamic range of 4-7 orders of magnitude (a bright
source core can be a million times the surface brightness of the disk halo).
Linearly mapping this to 256 display levels collapses everything to black or
white. The whole reason early plates look like beautiful nebulae and your
first matplotlib attempt looks like a black square is the stretch function.

WHAT A "STRETCH" IS
===================
A nonlinear remapping from data value to display intensity. The standard
choices, all in `astropy.visualization`:

  LinearStretch    — identity. Use only if data is already log-flux or similar.
  LogStretch       — emphasizes faint structure. Good for galaxy halos, jets.
  SqrtStretch      — middle ground. Standard for Hubble press images.
  AsinhStretch     — like log but well-behaved at zero. Standard for SDSS.
  PowerStretch     — gamma correction. Tunable.
  HistEqStretch    — histogram equalization. Maximum contrast, distorts levels.

Paired with an `Interval` (which sets the data range to map): MinMaxInterval,
PercentileInterval(99), ZScaleInterval (the SAOImage DS9 default, often best).

REFERENCES
==========
- Astropy stretches guide:
  https://docs.astropy.org/en/stable/visualization/normalization.html
- Lupton, Blanton, Fekete et al. 2004, "Preparing Red-Green-Blue Images"
  arXiv:astro-ph/0312483 — the seminal paper on asinh for color images.
- DS9 ZScale algorithm:
  http://iraf.net/forum/viewtopic.php?showtopic=139256
"""

from __future__ import annotations

from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import (
    AsinhStretch,
    ImageNormalize,
    LinearStretch,
    LogStretch,
    PercentileInterval,
    SqrtStretch,
    ZScaleInterval,
)
from matplotlib.figure import Figure

from .io import ImageData

StretchName = Literal["linear", "sqrt", "log", "asinh", "zscale"]


def _build_norm(array: np.ndarray, stretch: StretchName) -> ImageNormalize:
    """Build an ImageNormalize matching the named stretch.

    `zscale` is special: it's both an interval AND implies linear stretch
    (the convention from DS9). The others use percentile clipping (99.5%)
    to ignore extreme outlier pixels (cosmic rays, hot pixels) when setting
    the display range.
    """
    if stretch == "zscale":
        return ImageNormalize(array, interval=ZScaleInterval(),
                              stretch=LinearStretch())

    interval = PercentileInterval(99.5)
    stretches = {
        "linear": LinearStretch(),
        "sqrt": SqrtStretch(),
        "log": LogStretch(),
        # AsinhStretch's parameter `a` controls the linear-to-log transition.
        # a=0.1 is a reasonable default; smaller emphasizes faint stuff harder.
        "asinh": AsinhStretch(a=0.1),
    }
    return ImageNormalize(array, interval=interval, stretch=stretches[stretch])


def render_image(
    image: ImageData,
    *,
    stretch: StretchName = "asinh",
    cmap: str = "inferno",
    title: str | None = None,
    show_grid: bool = True,
    show_colorbar: bool = True,
    figsize: tuple[float, float] = (8, 8),
) -> Figure:
    """Render a single astronomical image with proper conventions.

    KEY CONVENTIONS APPLIED:
      - WCS axes if available: tick labels in RA/Dec sexagesimal, North up.
      - Asinh stretch by default (Lupton-style; works for almost everything).
      - 'inferno' colormap (perceptually uniform; replaces the old 'jet'
        which had artificial bands and is widely deprecated for science use).
      - Origin lower-left (FITS convention; numpy default is upper-left).

    PITFALL: If you forget origin='lower' your image is flipped vertically
    relative to every other tool (DS9, fv, aplpy, etc.) — sky is upside down.
    """
    fig = plt.figure(figsize=figsize, facecolor="#0e1117")

    if image.wcs is not None:
        ax = fig.add_subplot(111, projection=image.wcs)
        # Axis labels via WCSAxes — the projection knows what RA/Dec are.
        ax.set_xlabel("Right Ascension (J2000)", color="white")
        ax.set_ylabel("Declination (J2000)", color="white")
        if show_grid:
            ax.grid(color="white", ls=":", alpha=0.4)
        # Tick label color
        for axis in ("ra", "dec"):
            try:
                ax.coords[axis].set_ticklabel(color="white")
                ax.coords[axis].set_axislabel(axis.upper(), color="white")
            except Exception:
                pass
    else:
        ax = fig.add_subplot(111)
        ax.set_xlabel("X (pixels)", color="white")
        ax.set_ylabel("Y (pixels)", color="white")

    norm = _build_norm(image.array, stretch)
    im = ax.imshow(image.array, origin="lower", cmap=cmap, norm=norm,
                   interpolation="nearest")

    if show_colorbar:
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Counts" if image.array.dtype.kind in "ui" else "Flux (a.u.)",
                       color="white")
        cbar.ax.yaxis.set_tick_params(color="white")
        if cbar.outline is not None:
            cbar.outline.set_edgecolor("white")  # type: ignore[operator]
        for label in cbar.ax.get_yticklabels():
            label.set_color("white")

    if title:
        ax.set_title(title, color="white", fontsize=13, pad=12)

    # Dark theme to match Streamlit dark mode
    ax.set_facecolor("#0e1117")
    for spine in ax.spines.values():
        spine.set_color("white")

    fig.tight_layout()
    return fig


def render_event_image(
    binned_image: np.ndarray,
    extent: tuple[float, float, float, float],
    *,
    stretch: StretchName = "asinh",
    cmap: str = "inferno",
    title: str | None = None,
    energy_band_label: str | None = None,
    figsize: tuple[float, float] = (8, 8),
) -> Figure:
    """Render a binned X-ray event image (no WCS — sky-coord X/Y).

    Event-list images have sky-coordinate X,Y axes in detector units. WCS
    can in principle be constructed from event-file headers, but for v1
    we plot in those native units with a clear "X-ray sky pixels" label.
    """
    fig, ax = plt.subplots(figsize=figsize, facecolor="#0e1117")
    norm = _build_norm(binned_image, stretch)
    im = ax.imshow(binned_image, origin="lower", cmap=cmap, norm=norm,
                   extent=extent, interpolation="nearest")
    ax.set_xlabel("Sky X (detector pixels)", color="white")
    ax.set_ylabel("Sky Y (detector pixels)", color="white")
    if title:
        ax.set_title(title, color="white", fontsize=13, pad=12)
    if energy_band_label:
        ax.text(0.02, 0.98, energy_band_label, transform=ax.transAxes,
                color="white", fontsize=11, va="top",
                bbox=dict(facecolor="black", alpha=0.5, edgecolor="white", lw=0.5))

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Photon counts per pixel", color="white")
    cbar.ax.yaxis.set_tick_params(color="white")
    if cbar.outline is not None:
        cbar.outline.set_edgecolor("white")  # type: ignore[operator]
    for label in cbar.ax.get_yticklabels():
        label.set_color("white")

    ax.set_facecolor("#0e1117")
    for spine in ax.spines.values():
        spine.set_color("white")
    fig.tight_layout()
    return fig
