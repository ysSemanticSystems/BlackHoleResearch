"""
blackhole.provenance — reproducibility metadata attached to every figure.

Every ``render_*`` function in this package attaches a
``Provenance`` instance to ``fig.metadata["provenance"]``. The UI's
"Provenance" expander reads from it. ``save_figure(fig, path)`` writes
the PNG *plus* a sidecar JSON so an off-machine reader can reconstruct
which file, calibration version, and library version produced the plot.

Schema
------
fits_sha256       hex digest of the underlying FITS file (or "synthetic")
fits_path         resolved path of the source file
calibration_version    blackhole.calibration.CALIBRATION_VERSION at render
function_chain    tuple of function names that contributed
                  e.g. ("calibrate", "aperture_photometry_on", "render_sed")
library_version   blackhole.__version__
timestamp_utc     ISO-8601 UTC timestamp of the render call

References
----------
- ADR-005 (PHASE2_PLAN.md §6) — provenance contract.
- W3C PROV-O — long-term direction.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field, replace
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import matplotlib.figure

from . import __version__ as _lib_version
from .calibration import CALIBRATION_VERSION

SCHEMA_VERSION = "1"


@dataclass(frozen=True)
class Provenance:
    """Metadata bundle attached to every blackhole.* render output."""

    fits_sha256: str
    fits_path: str
    calibration_version: str
    function_chain: tuple[str, ...]
    library_version: str
    timestamp_utc: str
    schema_version: str = SCHEMA_VERSION
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a plain dict (sidecar JSON format)."""
        d = asdict(self)
        # Convert the tuple to a list so the JSON round-trip is symmetric.
        d["function_chain"] = list(self.function_chain)
        return d

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> Provenance:
        """Reconstruct a Provenance from a sidecar JSON dict."""
        return cls(
            fits_sha256=str(d["fits_sha256"]),
            fits_path=str(d["fits_path"]),
            calibration_version=str(d["calibration_version"]),
            function_chain=tuple(d.get("function_chain", ())),
            library_version=str(d["library_version"]),
            timestamp_utc=str(d["timestamp_utc"]),
            schema_version=str(d.get("schema_version", SCHEMA_VERSION)),
            extra=dict(d.get("extra", {})),
        )


# ---------------------------------------------------------------------------
# Construction helpers
# ---------------------------------------------------------------------------


def sha256_of_file(path: str | Path | None, chunk_size: int = 1 << 20) -> str:
    """Hex SHA-256 of a file's contents.

    Returns ``"synthetic"`` when ``path`` is None or unreadable — this is
    the convention for synthetic / in-memory inputs (tests, derived figures).
    """
    if path is None:
        return "synthetic"
    p = Path(path)
    if not p.exists() or not p.is_file():
        return "synthetic"
    h = hashlib.sha256()
    with p.open("rb") as f:
        while chunk := f.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()


def build_provenance(
    fits_path: str | Path | None,
    *,
    function_chain: tuple[str, ...],
    extra: dict[str, Any] | None = None,
) -> Provenance:
    """Construct a Provenance for the current render call.

    Parameters
    ----------
    fits_path
        Path to the FITS file (or None for purely synthetic input).
    function_chain
        Ordered tuple of function names that produced this figure, from
        the data load through the renderer. By convention the last entry
        is the renderer name (``"render_sed"`` etc.).
    extra
        Free-form key/value bag captured alongside the strict schema.
        Useful for stretch/cmap, aperture radii, fit ranges.
    """
    return Provenance(
        fits_sha256=sha256_of_file(fits_path),
        fits_path=str(Path(fits_path).resolve()) if fits_path else "synthetic",
        calibration_version=CALIBRATION_VERSION,
        function_chain=function_chain,
        library_version=_lib_version,
        timestamp_utc=datetime.now(UTC).isoformat(timespec="seconds"),
        extra=dict(extra or {}),
    )


def attach(fig: matplotlib.figure.Figure, provenance: Provenance) -> None:
    """Attach the Provenance to a Figure under ``fig.metadata['provenance']``.

    matplotlib's Figure doesn't have a public metadata dict, so we set
    one as an attribute (``fig.metadata``) if it isn't already present.
    """
    metadata: dict[str, Any] = getattr(fig, "metadata", None) or {}
    metadata["provenance"] = provenance
    fig.metadata = metadata  # type: ignore[attr-defined]


def get(fig: matplotlib.figure.Figure) -> Provenance | None:
    """Return the Provenance attached to a Figure, if any."""
    metadata = getattr(fig, "metadata", None)
    if not metadata:
        return None
    p = metadata.get("provenance")
    return p if isinstance(p, Provenance) else None


def extend_chain(provenance: Provenance, *names: str) -> Provenance:
    """Return a new Provenance with extra function names appended."""
    return replace(provenance, function_chain=provenance.function_chain + tuple(names))


# ---------------------------------------------------------------------------
# Save / load
# ---------------------------------------------------------------------------


def save_figure(
    fig: matplotlib.figure.Figure,
    path: str | Path,
    *,
    dpi: int = 150,
) -> tuple[Path, Path]:
    """Save the figure as PNG + sidecar JSON.

    Parameters
    ----------
    fig
        Figure with a Provenance attached via `attach`.
    path
        Output path. The PNG is written here (suffix is .png if absent).
        The sidecar JSON is written alongside with the same stem + .json.
    dpi
        PNG resolution.

    Returns
    -------
    (png_path, json_path)
        Filesystem paths of the two artefacts.

    Raises
    ------
    ValueError
        If the figure has no Provenance attached.
    """
    provenance = get(fig)
    if provenance is None:
        raise ValueError(
            "Figure has no Provenance attached. Call provenance.attach(fig, ...)"
            " in the renderer before saving."
        )
    p = Path(path)
    if p.suffix.lower() != ".png":
        p = p.with_suffix(".png")
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=dpi, facecolor=fig.get_facecolor())
    json_path = p.with_suffix(".json")
    json_path.write_text(json.dumps(provenance.to_dict(), indent=2))
    return p, json_path


def load_sidecar(json_path: str | Path) -> Provenance:
    """Load a Provenance back from the sidecar JSON."""
    text = Path(json_path).read_text()
    return Provenance.from_dict(json.loads(text))


# ---------------------------------------------------------------------------
# Human-readable export
# ---------------------------------------------------------------------------


def as_table_rows(provenance: Provenance) -> list[tuple[str, str]]:
    """Render the Provenance as a list of ``(field, value)`` strings for
    a Streamlit table or markdown dump."""
    rows: list[tuple[str, str]] = [
        ("fits_path",            provenance.fits_path),
        ("fits_sha256",          provenance.fits_sha256[:16] + "…"),
        ("function_chain",       " → ".join(provenance.function_chain) or "—"),
        ("library_version",      provenance.library_version),
        ("calibration_version",  provenance.calibration_version),
        ("timestamp_utc",        provenance.timestamp_utc),
        ("schema_version",       provenance.schema_version),
    ]
    for k, v in provenance.extra.items():
        rows.append((f"extra.{k}", str(v)))
    return rows


def as_bibtex_note(provenance: Provenance, target: str | None = None) -> str:
    """Compact citation note suitable for a clipboard 'copy citation' button."""
    pieces = [
        f"Figure generated by BlackHoleResearch v{provenance.library_version}",
        f"(calibration v{provenance.calibration_version})",
        f"on {provenance.timestamp_utc}",
        f"from {Path(provenance.fits_path).name}",
        f"(sha256={provenance.fits_sha256[:16]}…)",
    ]
    if target:
        pieces.insert(0, f"Target: {target}.")
    return " ".join(pieces) + "."
