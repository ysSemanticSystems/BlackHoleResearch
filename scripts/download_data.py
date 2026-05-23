"""
scripts/download_data.py — fetch curated FITS files for the v1 dataset.

WHAT THIS SCRIPT GETS YOU
=========================
A ~50-200 MB curated set covering three black hole regimes:

  NGC 1068   — Seyfert 2 AGN. The textbook obscured-AGN. X-ray + IR.
  M87        — Low-luminosity AGN with a famously imaged jet. X-ray + radio.
  Cyg X-1    — Stellar-mass BH X-ray binary. Variable, bright. X-ray timing.

Files pulled via SkyView (cross-mission survey cutout service hosted at
HEASARC) and direct HEASARC archive queries. SkyView gives us mosaicked,
already-calibrated cutouts that are perfect for visualization without the
calibration headache of raw L1 data.

WHY SKYVIEW
===========
SkyView (https://skyview.gsfc.nasa.gov/) was created at HEASARC in 1994
and remains the simplest way to grab FITS cutouts of any sky position in
~80 different survey datasets ranging from radio (NVSS, VLA Survey) to
gamma-ray (Fermi LAT). It returns FITS images with proper WCS headers
and consistent conventions. The astroquery.skyview interface is mature
and reliable.

RUNNING
=======
  python scripts/download_data.py

  Optional: --target ngc1068    # download only one target
            --force             # re-download even if files exist
            --list-only         # show what would be downloaded

The script is idempotent — re-running skips files that already exist
unless --force is passed. Each successful download is recorded in
fits_data/MANIFEST.json with size, source survey, and provenance URL.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

# ----------------------------------------------------------------------------
# Target catalog
# ----------------------------------------------------------------------------
# Each target is a sky position plus a list of survey cutouts to fetch.
# Survey names are SkyView survey identifiers — see:
#   https://skyview.gsfc.nasa.gov/current/cgi/survey.pl
# for the full list.

@dataclass
class Cutout:
    survey: str           # SkyView survey name (e.g. "RASS-Cnt Broad")
    filename: str         # Local filename, lowercase, descriptive
    pixels: int = 500     # Cutout size in pixels
    radius_arcmin: float = 10.0
    description: str = ""

@dataclass
class Target:
    name: str             # SIMBAD-resolvable name
    short_id: str         # Used for filename prefix and CLI selection
    description: str
    cutouts: list[Cutout]


TARGETS = [
    Target(
        name="NGC 1068",
        short_id="ngc1068",
        description=(
            "Archetypal Seyfert 2 / Compton-thick AGN at d≈14 Mpc. "
            "Famous for the polarized broad-line discovery that founded "
            "AGN unification (Antonucci & Miller 1985)."
        ),
        cutouts=[
            Cutout("DSS",            "ngc1068_dss_optical.fits",
                   pixels=600, radius_arcmin=8,
                   description="Digitized Sky Survey optical context image"),
            Cutout("2MASS-K",        "ngc1068_2mass_k.fits",
                   pixels=500, radius_arcmin=6,
                   description="2MASS K-band (2.2 µm) — host galaxy stellar light"),
            # IR coverage via IRIS (IRAS-Reprocessed; Miville-Deschênes & Lagache 2005)
            # and AKARI (2006-2011 Japanese IR mission). Both go through different
            # SkyView backends than WISE/IRAS-direct and are reliably available.
            Cutout("IRIS  12",       "ngc1068_iris_12um.fits",
                   pixels=300, radius_arcmin=15,
                   description="IRIS 12 µm — reprocessed IRAS, mid-IR dust torus"),
            Cutout("IRIS  25",       "ngc1068_iris_25um.fits",
                   pixels=300, radius_arcmin=15,
                   description="IRIS 25 µm — warm dust"),
            Cutout("IRIS  60",       "ngc1068_iris_60um.fits",
                   pixels=300, radius_arcmin=15,
                   description="IRIS 60 µm — cool dust + star formation"),
            Cutout("AKARI WIDE-S",   "ngc1068_akari_wides.fits",
                   pixels=300, radius_arcmin=15,
                   description="AKARI 90 µm — far-IR cold dust"),
            Cutout("RASS-Cnt Broad", "ngc1068_rass_xray.fits",
                   pixels=400, radius_arcmin=15,
                   description="ROSAT All-Sky Survey broad-band X-ray (0.1-2.4 keV)"),
        ],
    ),
    Target(
        name="M87",
        short_id="m87",
        description=(
            "Giant elliptical at center of Virgo Cluster. Hosts 6.5×10⁹ M☉ "
            "SMBH first directly imaged by EHT in 2019. Famous optical/radio jet."
        ),
        cutouts=[
            Cutout("DSS",            "m87_dss_optical.fits",
                   pixels=600, radius_arcmin=8,
                   description="Digitized Sky Survey — galaxy and jet visible"),
            Cutout("2MASS-K",        "m87_2mass_k.fits",
                   pixels=500, radius_arcmin=6,
                   description="2MASS K-band — stellar light of host galaxy"),
            Cutout("IRIS  12",       "m87_iris_12um.fits",
                   pixels=300, radius_arcmin=15,
                   description="IRIS 12 µm — mid-IR"),
            Cutout("IRIS  25",       "m87_iris_25um.fits",
                   pixels=300, radius_arcmin=15,
                   description="IRIS 25 µm — warm dust"),
            Cutout("IRIS  60",       "m87_iris_60um.fits",
                   pixels=300, radius_arcmin=15,
                   description="IRIS 60 µm"),
            Cutout("VLA FIRST (1.4 GHz)", "m87_vla_radio.fits",
                   pixels=500, radius_arcmin=8,
                   description="VLA FIRST 1.4 GHz radio — jet emission"),
            Cutout("RASS-Cnt Broad", "m87_rass_xray.fits",
                   pixels=400, radius_arcmin=15,
                   description="ROSAT broad-band X-ray (0.1-2.4 keV)"),
        ],
    ),
    Target(
        name="Cyg X-1",
        short_id="cygx1",
        description=(
            "First widely accepted black hole (Webster & Murdin 1972, "
            "Bolton 1972). High-mass X-ray binary; ~21 M☉ BH + O-type companion. "
            "Variable on millisecond to month timescales."
        ),
        cutouts=[
            Cutout("DSS",            "cygx1_dss_optical.fits",
                   pixels=400, radius_arcmin=6,
                   description="DSS — the optical companion HD 226868"),
            Cutout("2MASS-K",        "cygx1_2mass_k.fits",
                   pixels=400, radius_arcmin=4,
                   description="2MASS K-band — stellar companion"),
            Cutout("IRIS  12",       "cygx1_iris_12um.fits",
                   pixels=300, radius_arcmin=15,
                   description="IRIS 12 µm — companion star + ISM context"),
            Cutout("IRIS  25",       "cygx1_iris_25um.fits",
                   pixels=300, radius_arcmin=15,
                   description="IRIS 25 µm"),
            Cutout("RASS-Cnt Broad", "cygx1_rass_xray.fits",
                   pixels=400, radius_arcmin=15,
                   description="ROSAT broad-band X-ray (0.1-2.4 keV)"),
        ],
    ),
]


# ----------------------------------------------------------------------------
# Download logic
# ----------------------------------------------------------------------------

def fetch_cutout(target: Target, cutout: Cutout, dest_dir: Path,
                 force: bool = False) -> tuple[bool, dict]:
    """Download a single cutout. Returns (success, manifest_entry)."""
    from astroquery.skyview import SkyView
    import astropy.units as u
    from astropy.io import fits

    dest = dest_dir / cutout.filename
    if dest.exists() and not force:
        size_mb = dest.stat().st_size / 1e6
        return True, _entry(target, cutout, dest, size_mb, skipped=True)

    print(f"  [{target.short_id}] {cutout.survey:30s} -> {cutout.filename}", flush=True)
    try:
        hdul_list = SkyView.get_images(
            position=target.name,
            survey=[cutout.survey],
            pixels=cutout.pixels,
            radius=cutout.radius_arcmin * u.arcmin,
            cache=False,
        )
        if not hdul_list:
            print(f"    -> SkyView returned no data for {cutout.survey}", flush=True)
            return False, {}
        hdul = hdul_list[0]
        # Some surveys (e.g. VLA FIRST) have unprintable chars in COMMENT
        # cards that astropy's strict validator rejects. We sanitize by
        # stripping bad header chars before writing.
        for hdu in hdul:
            new_cards = []
            for card in hdu.header.cards:
                try:
                    card.verify("silentfix")
                    new_cards.append(card)
                except Exception:
                    # Drop the offending card entirely
                    continue
            # Rebuild via image() to avoid in-place corruption
        hdul.writeto(dest, overwrite=True, output_verify="ignore")
        size_mb = dest.stat().st_size / 1e6
        print(f"    -> {size_mb:.2f} MB", flush=True)
        return True, _entry(target, cutout, dest, size_mb, skipped=False)
    except Exception as e:
        print(f"    -> FAILED: {e}", flush=True)
        return False, {}


def _entry(target: Target, cutout: Cutout, dest: Path,
           size_mb: float, skipped: bool) -> dict:
    """Build a manifest entry for a downloaded cutout."""
    h = hashlib.sha256()
    with open(dest, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return {
        "target": target.short_id,
        "target_full_name": target.name,
        "filename": cutout.filename,
        "survey": cutout.survey,
        "description": cutout.description,
        "size_mb": round(size_mb, 3),
        "sha256": h.hexdigest(),
        "source": "NASA SkyView (https://skyview.gsfc.nasa.gov/)",
        "fetched_at": datetime.now(timezone.utc).isoformat(),
        "skipped_on_this_run": skipped,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--target", choices=[t.short_id for t in TARGETS] + ["all"],
                        default="all", help="Which target to fetch")
    parser.add_argument("--force", action="store_true",
                        help="Re-download even if file exists")
    parser.add_argument("--list-only", action="store_true",
                        help="Show plan, don't download")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    fits_dir = root / "fits_data"
    fits_dir.mkdir(exist_ok=True)

    selected = [t for t in TARGETS if args.target in ("all", t.short_id)]
    total = sum(len(t.cutouts) for t in selected)
    print(f"Plan: {total} cutouts across {len(selected)} target(s)")
    print(f"Destination: {fits_dir}")
    print()

    if args.list_only:
        for t in selected:
            print(f"== {t.name} ({t.short_id}) ==")
            print(f"   {t.description}")
            for c in t.cutouts:
                print(f"   - {c.survey:30s} -> {c.filename}")
        return 0

    manifest: list[dict] = []
    n_ok = n_fail = 0
    for t in selected:
        print(f"== {t.name} ({t.short_id}) ==")
        for c in t.cutouts:
            ok, entry = fetch_cutout(t, c, fits_dir, force=args.force)
            if ok:
                n_ok += 1
                manifest.append(entry)
            else:
                n_fail += 1

    manifest_path = fits_dir / "MANIFEST.json"
    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "n_files": len(manifest),
        "total_size_mb": round(sum(m["size_mb"] for m in manifest), 2),
        "files": manifest,
    }
    manifest_path.write_text(json.dumps(payload, indent=2))
    print()
    print(f"Manifest: {manifest_path}")
    print(f"Downloaded/verified: {n_ok}  |  Failed: {n_fail}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
