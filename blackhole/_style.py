"""
blackhole._style — single source of truth for matplotlib dark-mode styling.

Every renderer in `blackhole.*` is supposed to call ``apply_dark(ax)``
on each Axes it creates. Doing so:

- Keeps the white-on-dark colour palette uniform across plots.
- Makes future style changes a one-file edit rather than a grep-and-replace.
- Lets a single grep test enforce that we don't re-introduce the inline
  ``for spine in ax.spines.values(): spine.set_color("white")`` pattern
  in more than one module (cf. PHASE2_PLAN.md M6 exit criterion).
"""

from __future__ import annotations

import matplotlib.axes
import matplotlib.figure

DARK_BG = "#0e1117"        # matches Streamlit's default dark theme
DARK_FG = "white"
DARK_GRID = "white"
DARK_GRID_ALPHA = 0.3
DARK_GRID_LS = ":"


def apply_dark(ax: matplotlib.axes.Axes) -> None:
    """Apply the project's dark-mode styling to a single Axes.

    Idempotent. Safe to call multiple times on the same Axes.
    """
    ax.set_facecolor(DARK_BG)
    for spine in ax.spines.values():
        spine.set_color(DARK_FG)
    ax.tick_params(colors=DARK_FG)
    for label in (ax.xaxis.label, ax.yaxis.label, ax.title):
        label.set_color(DARK_FG)
    ax.grid(True, ls=DARK_GRID_LS, alpha=DARK_GRID_ALPHA, color=DARK_GRID)


def apply_dark_figure(fig: matplotlib.figure.Figure) -> None:
    """Apply dark background to a Figure (use after fig creation)."""
    fig.patch.set_facecolor(DARK_BG)


__all__ = [
    "DARK_BG", "DARK_FG", "DARK_GRID", "DARK_GRID_ALPHA", "DARK_GRID_LS",
    "apply_dark", "apply_dark_figure",
]
