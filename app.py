"""
app.py — Streamlit UI for BlackHoleResearch.

Launch:
    streamlit run app.py

Layout: compact top banner with a file selector and metric row, then tabs.
The Overview tab is the landing surface (thumbnail grid of available files).
Display-side controls (stretch, colormap) live inside the relevant tabs so
they sit next to the figures they configure.

The "science banner" deep-dive is deferred to Phase 2 M1 (source-catalog
work), at which point we can drive it from the active target rather than
a static expander block. Tracked in PHASE2_PLAN.md M1.
"""

from __future__ import annotations

import io as _stdio
import json
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import streamlit as st

from blackhole import calibration as calmod
from blackhole import catalog as cat
from blackhole import io as bhio
from blackhole import lightcurves as lc
from blackhole import photometry as phot
from blackhole import sed as sedmod
from blackhole import spectra as sp
from blackhole import wcs_plot as wp
from docs import literature_seds as litseds

# ---------------------------------------------------------------------------
# Page setup
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="BlackHoleResearch — Multiwavelength FITS Explorer",
    page_icon="🕳️",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
    .block-container { padding-top: 1rem; padding-bottom: 1rem; }
    h1 { margin-bottom: 0.25rem; }
    h1 + p { margin-top: 0; }
    .metric-card {
        background-color: #1a1d29; border: 1px solid #303446;
        border-radius: 8px; padding: 12px; margin-bottom: 10px;
    }
    section[data-testid="stSidebar"] { width: 280px !important; }
    </style>
    """,
    unsafe_allow_html=True,
)


# ---------------------------------------------------------------------------
# Compact header
# ---------------------------------------------------------------------------

st.title("BlackHoleResearch")
st.caption(
    "Multiwavelength FITS explorer for black hole astrophysics — "
    "NGC 1068, M87, Cyg X-1 across radio, IR, optical, and X-ray."
)


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

DATA_DIR = Path(__file__).resolve().parent / "fits_data"
MANIFEST_PATH = DATA_DIR / "MANIFEST.json"


def list_fits_files() -> list[Path]:
    return sorted(DATA_DIR.glob("*.fits"))


files = list_fits_files()
if not files:
    st.error(
        "No FITS files in `fits_data/`.\n\n"
        "From the project root:\n"
        "```bash\npython scripts/download_data.py\n```"
    )
    st.stop()

file_names = [f.name for f in files]


# Sidebar now carries only orientation material — file selection moved
# in-page so it stays adjacent to the figure it drives.
with st.sidebar:
    st.header("About this dataset")

    if MANIFEST_PATH.exists():
        manifest = json.loads(MANIFEST_PATH.read_text())
        st.caption(
            f"**{manifest.get('n_files', 0)}** files · "
            f"**{manifest.get('total_size_mb', 0)} MB**"
        )
        st.caption(f"Generated: `{manifest.get('generated_at', 'unknown')[:19]}`")

    st.markdown(
        "**Targets in this build**\n"
        "- **NGC 1068** — Compton-thick obscured Seyfert 2\n"
        "- **M87** — low-luminosity AGN with jet\n"
        "- **Cyg X-1** — stellar-mass X-ray binary\n\n"
        "Phase 2 will add **JWST Little Red Dots**.\n"
    )

    st.markdown("---")
    st.markdown(
        "**Data sources**\n"
        "- [NASA SkyView](https://skyview.gsfc.nasa.gov/)\n"
        "- [HEASARC](https://heasarc.gsfc.nasa.gov/)\n"
        "- [MAST](https://mast.stsci.edu/)\n"
    )

    st.markdown("---")
    st.markdown(
        "**Project**\n"
        "- [Repository](https://github.com/ysSemanticSystems/BlackHoleResearch)\n"
        "- [Phase 2 plan](https://github.com/ysSemanticSystems/BlackHoleResearch/blob/main/PHASE2_PLAN.md)\n"
        "- [Pitfalls catalog](https://github.com/ysSemanticSystems/BlackHoleResearch/blob/main/docs/PITFALLS.md)\n"
    )


# ---------------------------------------------------------------------------
# Helpers (defined before the page flow that calls them)
# ---------------------------------------------------------------------------

_TYPE_PRETTY: dict[str, str] = {
    "seyfert1":  "Seyfert 1 (unobscured)",
    "seyfert2":  "Seyfert 2 (obscured)",
    "llagn":     "Low-lum. AGN",
    "xrb_hmxb":  "HMXB",
    "xrb_lmxb":  "LMXB",
    "lrd":       "Little Red Dot",
    "quasar":    "Quasar",
    "blazar":    "Blazar",
    "tde":       "Tidal disruption",
}


def detect_source_from_filename(name: str) -> cat.Source | None:
    """Resolve the active FITS filename to a catalog Source, or None."""
    return cat.by_filename(name)


@st.cache_data(show_spinner=False)
def cached_inspect(path_str: str):
    return bhio.inspect(path_str)


@st.cache_data(show_spinner=False)
def thumbnail_png(path_str: str, stretch: str = "asinh", cmap: str = "inferno") -> bytes:
    """Render a small preview PNG for the Overview grid.

    Caches by (path, stretch, cmap) so file pickers stay fast. For files
    with no 2D image HDU we bin the EVENTS extension at low resolution.
    Returns PNG bytes ready for st.image; raises on totally invalid files
    so the caller can surface an inline error rather than blanking out
    the tab.
    """
    hdu_list = bhio.inspect(path_str)
    fig = None
    try:
        image_idx = next(
            (h.index for h in hdu_list
             if h.hdu_type in ("PrimaryHDU", "ImageHDU", "CompImageHDU")
             and h.shape and len(h.shape) == 2),
            None,
        )
        if image_idx is not None:
            img = bhio.load_image(path_str, hdu_index=image_idx)
            fig = wp.render_image(
                img, stretch=stretch, cmap=cmap, title=None,
                show_colorbar=False, figsize=(3.4, 3.4),
            )
        elif any(h.name.upper() == "EVENTS" for h in hdu_list):
            evlist = bhio.load_events(path_str)
            arr, extent = bhio.bin_to_image(evlist, bins=192)
            fig = wp.render_event_image(
                arr, extent, stretch=stretch, cmap=cmap,
                title=None, energy_band_label=None,
                figsize=(3.4, 3.4),
            )
        else:
            raise ValueError("no renderable HDU")

        buf = _stdio.BytesIO()
        fig.savefig(buf, format="png", dpi=110, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
        return buf.getvalue()
    finally:
        if fig is not None:
            plt.close(fig)


# ---------------------------------------------------------------------------
# In-page selector + metric strip
# ---------------------------------------------------------------------------


# Streamlit ordering note: selectbox state can be read directly via the
# default-return value, so this single widget drives the whole page.
sel_col, m1, m2, m3, m4 = st.columns([2.4, 1, 1, 1, 1])
selected_name = sel_col.selectbox(
    "FITS file",
    file_names,
    index=0,
    help="Files in fits_data/. The Overview tab below shows thumbnails for all of them.",
    label_visibility="visible",
)
selected_path = DATA_DIR / selected_name
hdus = cached_inspect(str(selected_path))
primary = hdus[0]

m1.metric("Telescope",  primary.telescope or "—")
m2.metric("Instrument", primary.instrument or "—")
m3.metric("HDUs",       len(hdus))
m4.metric("File size",  f"{selected_path.stat().st_size / 1e6:.1f} MB")

# Catalog-driven science banner. Renders only when the filename resolves
# to a catalogued source; otherwise we stay quiet so files outside the
# Phase-1 target set don't get bogus physics annotations.
active_source = detect_source_from_filename(selected_name)
if active_source is not None:
    L_edd = cat.eddington_luminosity_of(active_source).value
    d_quant = cat.distance_to(active_source)
    distance_str = (
        f"{d_quant.to(u.kpc).value:.2f} kpc"
        if d_quant.value < 0.1
        else f"{d_quant.value:.1f} Mpc"
    )

    if active_source.m_bh_msun is not None:
        if active_source.m_bh_msun >= 1e6:
            m_bh_str = f"{active_source.m_bh_msun / 1e6:.2f}×10⁶ M☉"
        else:
            m_bh_str = f"{active_source.m_bh_msun:.1f} M☉"
        if active_source.m_bh_msun >= 1e9:
            m_bh_str = f"{active_source.m_bh_msun / 1e9:.2f}×10⁹ M☉"
        if active_source.m_bh_err_msun:
            m_bh_str += f" ± {active_source.m_bh_err_msun:.2g}"
    else:
        m_bh_str = "—"

    src_cols = st.columns([1.4, 1, 1, 1, 1])
    src_cols[0].metric("Source", active_source.name)
    src_cols[1].metric("Type", _TYPE_PRETTY.get(active_source.type, active_source.type))
    src_cols[2].metric("M_BH", m_bh_str)
    src_cols[3].metric("Distance", distance_str)
    src_cols[4].metric("L_Edd (erg/s)", f"{L_edd:.2e}")

    with st.expander(f"References for {active_source.name}"):
        z_line = (
            f"- **Redshift**: z = {active_source.redshift:.5f}"
            if active_source.redshift is not None
            else "- **Redshift**: — (Galactic source)"
        )
        st.markdown(
            f"- **M_BH ref**: {active_source.m_bh_ref}\n"
            f"- **Distance ref**: {active_source.distance_ref}\n"
            f"{z_line}\n"
            f"- **Notes**: {active_source.notes or '—'}\n"
        )

st.divider()


# ---------------------------------------------------------------------------
# Main tabs
# ---------------------------------------------------------------------------

(
    tab_overview,
    tab_image,
    tab_sed,
    tab_xray,
    tab_spectrum,
    tab_lightcurve,
    tab_inspect,
    tab_about,
) = st.tabs([
    "Overview",
    "Image",
    "SED Builder",
    "X-ray Events",
    "Spectrum",
    "Light Curve",
    "Inspect",
    "About",
])

# Overview is the landing tab. Image / X-ray Events configure their own
# stretch+cmap in per-tab expanders because changing them mid-analysis
# should not require a sidebar round-trip.

# --------------------------- Overview tab ---------------------------------
with tab_overview:
    st.subheader("Dataset at a glance")
    st.caption(
        "One thumbnail per FITS file in `fits_data/`. Click the file name "
        "in the selector at the top of the page to drill into any of them."
    )

    overview_cmap = st.selectbox(
        "Thumbnail colormap",
        ["inferno", "viridis", "magma", "plasma", "gray"],
        index=0,
        key="overview_cmap",
        help="Affects the Overview grid only; per-tab views keep their own colormaps.",
    )

    n_cols = 3
    cols = st.columns(n_cols)
    for i, path in enumerate(files):
        col = cols[i % n_cols]
        with col:
            try:
                png = thumbnail_png(str(path), stretch="asinh", cmap=overview_cmap)
                st.image(png, caption=path.name, use_container_width=True)
            except Exception as e:
                st.warning(f"{path.name}: {e}")
            # Compact header line for at-a-glance context.
            try:
                preview_hdus = cached_inspect(str(path))
                p0 = preview_hdus[0]
                badges = " · ".join(filter(None, [
                    p0.telescope, p0.instrument,
                    f"{len(preview_hdus)} HDU{'s' if len(preview_hdus) != 1 else ''}",
                ]))
                if badges:
                    st.caption(badges)
            except Exception:
                pass

# --------------------------- Inspect tab ----------------------------------
with tab_inspect:
    st.subheader("HDU structure")
    st.caption(
        "Each row is one HDU. Image HDUs hold 2D arrays; BinTable HDUs hold "
        "structured rows (event lists, spectra, light curves). The Shape column "
        "shows row counts for tables and pixel dimensions for images."
    )
    import pandas as pd
    df = pd.DataFrame([
        {
            "Index": h.index,
            "Name": h.name,
            "Type": h.hdu_type,
            "Shape": str(h.shape) if h.shape else "—",
            "Columns": ", ".join(h.columns) if h.columns else "—",
            "Telescope": h.telescope or "—",
            "Instrument": h.instrument or "—",
            "Exposure (s)": f"{h.exposure_s:.1f}" if h.exposure_s else "—",
        }
        for h in hdus
    ])
    st.dataframe(df, use_container_width=True, hide_index=True)

    st.subheader("Primary header")
    st.caption("Selected keywords — full header has more.")
    hdr = bhio.header_dict(str(selected_path), 0)
    interesting = ["NAXIS", "NAXIS1", "NAXIS2", "BUNIT", "EXPTIME",
                   "TELESCOP", "INSTRUME", "OBJECT", "DATE-OBS", "DATE",
                   "EQUINOX", "RADECSYS", "RA", "DEC", "CTYPE1", "CTYPE2",
                   "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "BMAJ", "BMIN",
                   "ORIGIN", "AUTHOR", "REFERENC"]
    rows = [{"Keyword": k, "Value": str(hdr[k])} for k in interesting if k in hdr]
    if rows:
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
    else:
        st.info("None of the standard keywords were present in HDU 0. Check secondary HDUs.")

# --------------------------- Image tab ----------------------------------
with tab_image:
    st.subheader("Sky image with WCS")
    st.caption(
        "Renders 2D Image HDUs with the World Coordinate System overlaid. "
        "Axis labels become RA/Dec in J2000 sexagesimal when WCS is present."
    )

    with st.expander("Display options", expanded=False):
        oc1, oc2 = st.columns(2)
        image_stretch = oc1.selectbox(
            "Stretch",
            ["asinh", "log", "sqrt", "linear", "zscale"],
            index=0,
            key="image_stretch",
            help="Nonlinear stretch for high-dynamic-range data. "
                 "asinh is the modern default; zscale matches SAOImage DS9.",
        )
        image_cmap = oc2.selectbox(
            "Colormap",
            ["inferno", "viridis", "magma", "plasma", "gray", "cubehelix", "twilight"],
            index=0,
            key="image_cmap",
            help="Perceptually uniform colormaps.",
        )

    image_hdu_idx = next(
        (h.index for h in hdus
         if h.hdu_type in ("PrimaryHDU", "ImageHDU", "CompImageHDU")
         and h.shape and len(h.shape) == 2),
        None,
    )

    if image_hdu_idx is None:
        st.warning(
            "No 2D image HDU found in this file. If it's an X-ray event list, "
            "use the **X-ray Events** tab to bin it into an image first."
        )
    else:
        try:
            img = bhio.load_image(str(selected_path), hdu_index=image_hdu_idx)
            wcs_status = "WCS detected ✓" if img.wcs is not None else "no WCS — pixel coordinates only"
            st.caption(f"HDU {image_hdu_idx} · shape {img.array.shape} · {wcs_status}")
            fig = wp.render_image(
                img, stretch=image_stretch, cmap=image_cmap,
                title=f"{selected_name}  ·  {primary.telescope or ''} {primary.instrument or ''}".strip(),
            )
            st.pyplot(fig, use_container_width=True)

            with st.expander("Quick statistics"):
                import numpy as np
                arr = img.array
                st.write({
                    "min": float(np.nanmin(arr)),
                    "max": float(np.nanmax(arr)),
                    "mean": float(np.nanmean(arr)),
                    "median": float(np.nanmedian(arr)),
                    "std": float(np.nanstd(arr)),
                    "shape": arr.shape,
                    "dtype": str(arr.dtype),
                })
        except Exception as e:
            st.error(f"Could not render image: {e}")

# --------------------------- X-ray Events tab ----------------------------------
with tab_xray:
    st.subheader("X-ray event list → binned image")
    st.caption(
        "Bins photon event lists (Chandra / XMM-Newton / NuSTAR) into a 2D image. "
        "The energy filter separates hard (>2 keV) from soft (<2 keV) bands, which "
        "often look very different — hard X-rays penetrate obscuration that absorbs "
        "soft X-rays."
    )

    # Detect if this file has an EVENTS extension
    has_events = any(h.name.upper() == "EVENTS" for h in hdus)

    if not has_events:
        st.info(
            "This file doesn't appear to contain an X-ray event list "
            "(no EVENTS extension). Try a Chandra `evt2.fits`, "
            "XMM-Newton EPIC event file, or NuSTAR event list."
        )
    else:
        with st.expander("Display options", expanded=False):
            oc1, oc2 = st.columns(2)
            xray_stretch = oc1.selectbox(
                "Stretch",
                ["asinh", "log", "sqrt", "linear", "zscale"],
                index=0,
                key="xray_stretch",
            )
            xray_cmap = oc2.selectbox(
                "Colormap",
                ["inferno", "viridis", "magma", "plasma", "gray", "cubehelix"],
                index=0,
                key="xray_cmap",
            )

        c1, c2 = st.columns(2)
        bin_size = c1.slider("Bins per side", 128, 1024, 512, step=64)
        e_lo = c2.slider("Soft energy cutoff (eV)", 100, 3000, 500, step=100)
        e_hi = c2.slider("Hard energy cutoff (eV)", 1000, 12000, 7000, step=500)

        try:
            evlist = bhio.load_events(str(selected_path))
            st.caption(
                f"Loaded {evlist.times.size:,} photons · mission: {evlist.mission} · "
                f"energy unit: {evlist.energy_unit}"
            )
            image_arr, extent = bhio.bin_to_image(
                evlist, bins=bin_size,
                energy_range_ev=(e_lo, e_hi) if evlist.energy_unit == "eV" else None,
            )
            fig = wp.render_event_image(
                image_arr, extent,
                stretch=xray_stretch, cmap=xray_cmap,
                title=f"{selected_name}  ·  {evlist.mission}",
                energy_band_label=(f"{e_lo/1000:.1f}-{e_hi/1000:.1f} keV"
                                   if evlist.energy_unit == "eV" else None),
            )
            st.pyplot(fig, use_container_width=True)
        except Exception as e:
            st.error(f"Could not bin event image: {e}")

# --------------------------- Spectrum tab ----------------------------------
with tab_spectrum:
    st.subheader("1D X-ray spectrum")
    st.caption(
        "Plots an OGIP PHA spectrum and optionally fits a power law "
        "N(E) ∝ E^(−Γ). Without RMF/ARF response files the fit is in "
        "channel space — descriptive only, not calibrated."
    )
    with st.expander("Photon-index interpretation"):
        st.markdown(
            "- **Γ ≈ 1.7–2.1** → typical Type 1 (unobscured) AGN\n"
            "- **Γ < 1.7** → hard, reflection-dominated or obscured\n"
            "- **Γ > 2.5** → very soft; sometimes super-Eddington or tidal disruption event\n\n"
            "Production-grade fits require XSPEC or Sherpa with the appropriate "
            "response files."
        )

    has_spectrum = any(h.name.upper() == "SPECTRUM" for h in hdus)

    if not has_spectrum:
        st.info(
            "This file doesn't appear to contain a SPECTRUM HDU (OGIP PHA format). "
            "Try an X-ray spectrum file (often named `*_spec.pha` or `*.pi`)."
        )
    else:
        do_fit = st.checkbox("Fit a power law", value=True)
        try:
            spec = sp.load_pha_spectrum(str(selected_path))
            fit_result = None
            if do_fit:
                try:
                    g, ge, fc = sp.fit_power_law(spec)
                    fit_result = (g, ge, fc)
                    from blackhole.physics.spectral_xray import classify_photon_index
                    st.success(
                        f"**Γ = {g:.2f} ± {ge:.2f}** — *{classify_photon_index(g)}*"
                    )
                except Exception as fe:
                    st.warning(f"Fit failed: {fe}")
            fig = sp.render_spectrum(
                spec, fit=fit_result,
                title=f"{selected_name}  ·  {spec.mission} {spec.instrument}",
            )
            st.pyplot(fig, use_container_width=True)
        except Exception as e:
            st.error(f"Could not render spectrum: {e}")

# --------------------------- SED tab ----------------------------------
with tab_sed:
    st.subheader("Multi-band SED builder")
    st.caption(
        "Assembles a Spectral Energy Distribution across wavebands, with "
        "markers color-coded by band. Resolves the radio jet synchrotron "
        "component, the IR dust torus, the optical/UV disk thermal emission, "
        "and the X-ray corona power law on a single plot."
    )
    target_names = [s.name for s in cat.CATALOG]
    auto_source = detect_source_from_filename(selected_name)
    default_idx = (
        target_names.index(auto_source.name)
        if auto_source is not None and auto_source.name in target_names
        else 0
    )
    target_choice = st.selectbox(
        "Target",
        target_names,
        index=default_idx,
        help=(
            f"Auto-detected from filename → **{auto_source.name}**. Switch manually if needed."
            if auto_source is not None
            else "Filename doesn't match a catalogued target. Choose manually."
        ),
    )
    target_source = cat.by_name(target_choice)

    sed_controls_col1, sed_controls_col2 = st.columns([2, 1])
    with sed_controls_col1:
        use_local = st.checkbox(
            "Aperture-photometer the local FITS cutouts for this target",
            value=True,
            help=(
                "Runs `blackhole.photometry.aperture_photometry_on` against "
                "the cutouts in `fits_data/` whose filename mentions this "
                "target. Only files that calibrate (2MASS, IRIS, AKARI, "
                "FIRST, RASS) contribute; DSS raises and is skipped."
            ),
        )
        overlay_lit = st.checkbox(
            "Overlay published literature values for this target",
            value=True,
            help="Faint markers from `docs/literature_seds.py` for comparison.",
        )
    with sed_controls_col2:
        overlay_quasar = st.checkbox(
            "Overlay Elvis+1994 mean quasar template",
            value=False,
            help="Schematic broadband AGN spectrum for visual context only.",
        )

    sed_obj = sedmod.SED(target_name=target_choice)

    # ---------------- Local-cutout photometry (M2+M3 pipeline) -------------
    measured_rows: list[dict[str, object]] = []
    if use_local and target_source is not None:
        for path in sorted(list_fits_files(DATA_DIR)):
            ts = cat.by_filename(path.name)
            if ts is None or ts.name != target_source.name:
                continue
            try:
                img = bhio.load_image(path)
            except Exception as exc:
                measured_rows.append({"file": path.name, "status": f"load failed: {exc}"})
                continue
            try:
                cal_img = calmod.calibrate(img)
            except calmod.UncalibratedDataError as exc:
                measured_rows.append({"file": path.name, "status": f"skipped (uncalibrated): {exc}"})
                continue
            radius = phot.aperture_for_band(cal_img.band)
            try:
                res = phot.aperture_photometry_on(
                    cal_img,
                    target_source.coord,
                    aperture_radius=radius,
                    annulus_inner=radius * 1.5,
                    annulus_outer=radius * 2.5,
                    label=f"{cal_img.survey}",
                )
            except Exception as exc:
                measured_rows.append({"file": path.name, "status": f"photometry failed: {exc}"})
                continue

            # Wavelength tagging by survey identifier (so the SED renderer
            # can place the point in frequency space). Pure-best-effort
            # mapping; downstream Phase 2 will read this from a band registry.
            wavelength_um = {
                "2MASS-K":  2.159,
                "IRIS_12":  12.0,
                "IRIS_25":  25.0,
                "IRIS_60":  60.0,
                "IRIS_100": 100.0,
                "AKARI":     90.0,
                "VLA-FIRST": 214000.0,
                "RASS-broad": None,   # X-ray: handled via energy
            }.get(cal_img.survey)
            point = res.to_sed_point()
            if wavelength_um is not None:
                point.wavelength = wavelength_um * u.micron
            elif cal_img.survey == "RASS-broad":
                point.energy = 1.0 * u.keV  # representative broad-band energy
            sed_obj.add(point)
            measured_rows.append({
                "file": path.name,
                "survey": cal_img.survey,
                "flux": f"{res.flux:.3g}",
                "err":  f"{res.flux_err:.3g}",
                "upper_limit": res.upper_limit,
                "aperture_arcsec": float(res.aperture_radius.value),
            })

    # ---------------- Literature overlay (visual comparison only) ----------
    if overlay_lit:
        lit_sed, lit_xray = litseds.literature_for(target_choice)
        for label, wl_um, fnu_jy, band, src in lit_sed:
            sed_obj.add(sedmod.SEDPoint(
                label=f"[lit] {label}", wavelength=wl_um * u.micron,
                flux_density=fnu_jy * u.Jy, band=band,
                source=f"literature · {src}",
            ))
        for energy_keV, nfn, label, src in lit_xray:
            sed_obj.add(sedmod.SEDPoint(
                label=f"[lit] {label}", energy=energy_keV * u.keV,
                nu_f_nu=nfn * (u.erg / u.s / u.cm**2),
                band="xray", source=f"literature · {src}",
            ))

    if len(sed_obj.points) == 0:
        st.info(
            "No SED points to plot for this target. Enable the literature "
            "overlay or add calibratable cutouts (2MASS, IRIS, AKARI, FIRST, "
            "RASS) to `fits_data/` and re-run."
        )
    else:
        fig = sedmod.render_sed(sed_obj, overplot_quasar_template=overlay_quasar)
        st.pyplot(fig, use_container_width=True)

    with st.expander("Local photometry diagnostics"):
        if measured_rows:
            import pandas as pd
            st.dataframe(
                pd.DataFrame(measured_rows),
                use_container_width=True, hide_index=True,
            )
            st.caption(
                "Local aperture photometry from `blackhole.photometry`. "
                "Apertures size from `aperture_for_band(band)`; sky annulus "
                "is 1.5x → 2.5x the aperture. Upper limits are the 3σ "
                "background-noise threshold (see photutils docs)."
            )
        else:
            st.caption("No local cutouts contributed (try enabling the option above).")

    with st.expander("Literature points and references"):
        import pandas as pd
        lit_sed, lit_xray = litseds.literature_for(target_choice)
        rows = (
            [{"Label": lbl, "Band": band, "Source": src}
             for (lbl, _w, _f, band, src) in lit_sed]
            +
            [{"Label": lbl, "Band": "xray", "Source": src}
             for (_E, _nfn, lbl, src) in lit_xray]
        )
        if rows:
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
        else:
            st.caption("No literature points catalogued for this target.")
        st.caption(
            "Literature values are an opt-in overlay only. The data points "
            "you see on the plot when 'Aperture-photometer the local FITS "
            "cutouts' is enabled are computed locally."
        )

# --------------------------- Light Curve tab ----------------------------------
with tab_lightcurve:
    st.subheader("X-ray light curve and Lomb-Scargle periodogram")
    st.caption(
        "Bins photon arrival times into a uniform time-step light curve and "
        "computes a Lomb-Scargle periodogram for variability timescales. "
        "Useful for QPO searches and source-state characterization. Requires "
        "an event list (same input as the X-ray Events tab)."
    )

    has_events = any(h.name.upper() == "EVENTS" for h in hdus)
    if not has_events:
        st.info("No EVENTS HDU in this file. Use an X-ray event list.")
    else:
        bin_s = st.slider("Bin size (seconds)", 1, 1000, 100)
        try:
            evlist = bhio.load_events(str(selected_path))
            lc_obj = lc.bin_events_to_lightcurve(evlist, bin_size_s=bin_s)
            fig_lc = lc.render_lightcurve(lc_obj)
            st.pyplot(fig_lc, use_container_width=True)

            from blackhole.physics.variability import fractional_rms, fractional_rms_error
            fvar = fractional_rms(lc_obj.rates, lc_obj.errors)
            fvar_err = fractional_rms_error(lc_obj.rates, lc_obj.errors)
            if fvar == fvar:  # not NaN
                st.success(
                    f"**F_var = {fvar*100:.1f}% ± {fvar_err*100:.1f}%**  "
                    f"(fractional rms variability amplitude)"
                )
            else:
                st.info("F_var not well-defined (variability < noise level).")

            st.subheader("Lomb-Scargle periodogram")
            freqs, power = lc.lomb_scargle_periodogram(lc_obj)
            fig_p = lc.render_periodogram(freqs, power)
            st.pyplot(fig_p, use_container_width=True)
        except Exception as e:
            st.error(f"Light-curve computation failed: {e}")

# --------------------------- About tab ----------------------------------
with tab_about:
    st.subheader("About this project")
    st.markdown(
        """
### Phase 1 — Classical multiwavelength viewer (current build)
Three source regimes: NGC 1068 (obscured Seyfert 2), M87 (low-luminosity
AGN with jet), Cyg X-1 (stellar-mass black hole binary). Standard plots:
image with WCS, X-ray events, spectrum, SED, light curve.

### Phase 2 — Little Red Dots module
The JWST research frontier. Adds:
- LRD targets from JADES, CEERS, and COSMOS-Web (Greene et al. 2024, Maiolino et al. 2024)
- L_X / L_UV comparison against the Lusso & Risaliti 2016 relation that LRDs systematically violate
- NIRSpec spectral diagnostics (Balmer break detection, line ratio fits)
- IR-faintness diagnostic (Stern and Donley wedges with LRDs marked as outliers)

### Phase 3 — Calibrated spectral fitting
Sherpa or PyXspec integration with RMF/ARF handling and multi-component
models (power law + reflection + thermal disk).

### Data sources
| Mission / survey | Era | Use in this project |
|---|---|---|
| ROSAT All-Sky Survey | 1990–1991 | Broad-band X-ray cutouts via SkyView |
| 2MASS | 1997–2001 | Near-IR J/H/K cutouts |
| DSS | 1950s–1990s plates | Optical context images |
| WISE / AllWISE | 2010–2011 | Mid-IR (3.4 / 4.6 / 12 / 22 µm) |
| Chandra | 1999–present | X-ray events and spectra (via HEASARC) |
| NuSTAR | 2012–present | Hard X-ray (3–79 keV) |
| JWST | 2021–present | Phase 2 LRD photometry and NIRSpec |

### Standards
- **FITS Standard v4.0** (IAU FITS Working Group, 2018)
- **OGIP/92-007** for X-ray spectra
- **VOTable** and **WCS** for coordinates
- **astropy** 6.1+ for all I/O

### Pitfalls referenced in source
1. Event-file-as-image confusion → `io.py`
2. WCS missing or partial → `wcs_plot.py`
3. Linear stretch on high-dynamic-range data → `wcs_plot.py`
4. Unit confusion across bands (Jy / mag / keV / Hz) → `sed.py`
5. PHA spectrum requires RMF/ARF for a true energy axis → `spectra.py`
6. Big-endian byteswap and memmap invalidation → `io.py`
7. Light-curve gaps from unapplied GTI → `lightcurves.py`
8. WISE colors in Vega versus AB → `physics/infrared.py`

Full catalog in `docs/PITFALLS.md`.

### Source
[github.com/ysSemanticSystems/BlackHoleResearch](https://github.com/ysSemanticSystems/BlackHoleResearch)
        """
    )
