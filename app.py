"""
app.py — Streamlit UI for BlackHoleResearch.

Launch:
    streamlit run app.py

Layout: compact top banner with a metric row, then tabs for each
visualization. Background information is collected in collapsible
expanders so the visual content remains the focus.
"""

from __future__ import annotations

import json
from pathlib import Path

import streamlit as st

from blackhole import io as bhio
from blackhole import lightcurves as lc
from blackhole import sed as sedmod
from blackhole import spectra as sp
from blackhole import wcs_plot as wp


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

with st.expander("Background — what this tool does and why", expanded=False):
    st.markdown(
        """
        Black holes are detected through the radiation of accreting gas,
        which spans the electromagnetic spectrum. Each waveband traces a
        different physical region:

        - **Radio / sub-mm** — relativistic jets, synchrotron emission
        - **Infrared** — dusty torus around the central engine
        - **Optical / UV** — accretion disk (Shakura–Sunyaev thermal)
        - **X-ray** — hot corona, disk reflection, jet base
        - **Gamma-ray** — high-energy jet particles, inverse Compton

        This tool loads FITS data from multiple missions, displays each
        with appropriate stretches and WCS overlays, and assembles the
        **Spectral Energy Distribution** that ties the bands together.

        The current dataset covers three regimes: a Compton-thick
        obscured AGN (**NGC 1068**), a low-luminosity supermassive
        black hole with an extended jet (**M87**), and a stellar-mass
        black hole X-ray binary (**Cyg X-1**). Phase 2 will add
        **Little Red Dots** — the JWST-discovered population that
        appears infrared- and X-ray-weak despite hosting accreting
        black holes, in tension with the AGN unification model.

        **FITS format note:** A FITS file is a list of HDUs
        (Header/Data Units), each containing ASCII metadata and either
        an image array or a structured table. X-ray data from Chandra,
        XMM-Newton, and NuSTAR arrives as **event lists** (one row per
        detected photon) rather than 2D images. The Inspect tab shows
        the HDU structure for any loaded file; binning event lists into
        images is handled automatically.
        """
    )

st.divider()


# ---------------------------------------------------------------------------
# Sidebar — file selection and global controls
# ---------------------------------------------------------------------------

DATA_DIR = Path(__file__).resolve().parent / "fits_data"
MANIFEST_PATH = DATA_DIR / "MANIFEST.json"


def list_fits_files() -> list[Path]:
    return sorted(DATA_DIR.glob("*.fits"))


with st.sidebar:
    st.header("Controls")

    files = list_fits_files()
    if not files:
        st.error(
            "No FITS files in `fits_data/`.\n\n"
            "From the project root:\n"
            "```bash\npython scripts/download_data.py\n```"
        )
        st.stop()

    file_names = [f.name for f in files]
    selected_name = st.selectbox(
        "FITS file",
        file_names,
        help="Files in the fits_data/ directory.",
    )
    selected_path = DATA_DIR / selected_name

    st.divider()
    st.subheader("Image rendering")
    stretch = st.selectbox(
        "Stretch",
        ["asinh", "log", "sqrt", "linear", "zscale"],
        index=0,
        help="Nonlinear stretch for high-dynamic-range data. "
             "asinh is the modern default; zscale matches SAOImage DS9.",
    )
    cmap = st.selectbox(
        "Colormap",
        ["inferno", "viridis", "magma", "plasma", "gray", "cubehelix", "twilight"],
        index=0,
        help="Perceptually uniform colormaps.",
    )

    st.divider()
    if MANIFEST_PATH.exists():
        manifest = json.loads(MANIFEST_PATH.read_text())
        st.caption(
            f"Dataset: **{manifest.get('n_files', 0)}** files · "
            f"**{manifest.get('total_size_mb', 0)} MB**"
        )
        st.caption(f"Generated: `{manifest.get('generated_at', 'unknown')[:19]}`")

    st.caption(
        "Sources: [NASA SkyView](https://skyview.gsfc.nasa.gov/) · "
        "[HEASARC](https://heasarc.gsfc.nasa.gov/) · "
        "[MAST](https://mast.stsci.edu/)"
    )


# ---------------------------------------------------------------------------
# File metadata banner
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner=False)
def cached_inspect(path_str: str):
    return bhio.inspect(path_str)


hdus = cached_inspect(str(selected_path))
primary = hdus[0]

col_a, col_b, col_c, col_d = st.columns(4)
col_a.metric("Telescope", primary.telescope or "—")
col_b.metric("Instrument", primary.instrument or "—")
col_c.metric("HDUs", len(hdus))
col_d.metric("File size",
             f"{selected_path.stat().st_size / 1e6:.1f} MB")


# ---------------------------------------------------------------------------
# Main tabs
# ---------------------------------------------------------------------------

tab_image, tab_sed, tab_xray, tab_spectrum, tab_lightcurve, tab_inspect, tab_about = st.tabs([
    "Image",
    "SED Builder",
    "X-ray Events",
    "Spectrum",
    "Light Curve",
    "Inspect",
    "About",
])

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

    # Find first 2D image HDU
    image_hdu_idx = None
    for h in hdus:
        if h.hdu_type in ("PrimaryHDU", "ImageHDU", "CompImageHDU") and h.shape and len(h.shape) == 2:
            image_hdu_idx = h.index
            break

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
                img, stretch=stretch, cmap=cmap,
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
                stretch=stretch, cmap=cmap,
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
    with st.expander("Note on photometry source"):
        st.markdown(
            "Values are illustrative literature photometry, drawn from "
            "published catalog measurements (NED, AllWISE, 2MASS XSC, and "
            "the references shown in the data table). Aperture photometry "
            "extracted directly from the local FITS cutouts is planned for "
            "Phase 2."
        )

    target_choice = st.selectbox(
        "Target", ["NGC 1068", "M87", "Cyg X-1"], index=0,
    )

    # Hardcoded illustrative SED values per target. In a Phase 2 build,
    # these would be computed from the actual FITS cutouts via aperture
    # photometry. For now they're from published catalogs.
    #
    # Sources for each:
    #   NGC 1068: NED photometry compilation (https://ned.ipac.caltech.edu/)
    #   M87:      NED photometry; X-ray from Wilson & Yang 2002
    #   Cyg X-1:  Mid-IR from Mirabel et al. 1996; X-ray from Wilms+2006
    SEDS = {
        "NGC 1068": [
            # (label, wavelength_um, flux_density_jy, band, source)
            ("VLA 1.4 GHz",       214000.0,  2.0,    "radio",  "NED"),
            ("VLA 5 GHz",         60000.0,   1.4,    "radio",  "NED"),
            ("Spitzer 24 µm",     24.0,      20.0,   "ir",     "Bendo+2012"),
            ("WISE W4 22 µm",     22.0,      18.0,   "ir",     "AllWISE"),
            ("WISE W3 12 µm",     12.0,      18.0,   "ir",     "AllWISE"),
            ("Spitzer 8 µm",      8.0,       7.0,    "ir",     "Bendo+2012"),
            ("WISE W2 4.6 µm",    4.6,       2.6,    "ir",     "AllWISE"),
            ("WISE W1 3.4 µm",    3.4,       1.6,    "ir",     "AllWISE"),
            ("2MASS K 2.2 µm",    2.2,       1.0,    "ir",     "2MASS XSC"),
            ("2MASS J 1.25 µm",   1.25,      0.55,   "ir",     "2MASS XSC"),
            ("DSS R 0.66 µm",     0.66,      0.07,   "opt",    "DSS"),
            ("Swift UVOT UVW1",   0.26,      0.005,  "uv",     "Swift UVOT"),
            # X-ray points expressed in keV (handled differently below)
        ],
        "M87": [
            ("VLA 1.4 GHz core",  214000.0,  4.0,    "radio",  "NED"),
            ("VLA 5 GHz core",    60000.0,   2.9,    "radio",  "NED"),
            ("Spitzer 24 µm",     24.0,      0.4,    "ir",     "Shi+2007"),
            ("WISE W3 12 µm",     12.0,      0.5,    "ir",     "AllWISE"),
            ("WISE W1 3.4 µm",    3.4,       0.5,    "ir",     "AllWISE"),
            ("2MASS K 2.2 µm",    2.2,       1.5,    "ir",     "2MASS"),
            ("DSS R 0.66 µm",     0.66,      0.6,    "opt",    "DSS"),
        ],
        "Cyg X-1": [
            ("2MASS K 2.2 µm",    2.2,       0.8,    "ir",     "2MASS"),
            ("2MASS J 1.25 µm",   1.25,      0.55,   "ir",     "2MASS"),
            ("DSS R 0.66 µm",     0.66,      0.4,    "opt",    "DSS"),
        ],
    }

    # X-ray points are in keV
    XRAY = {
        "NGC 1068":  [(2.0, 5e-13, "Chandra 2 keV", "Bauer+2015"),
                      (10.0, 1e-12, "NuSTAR 10 keV", "Marinucci+2016"),
                      (30.0, 5e-12, "NuSTAR 30 keV", "Marinucci+2016")],
        "M87":       [(1.0, 5e-12, "Chandra 1 keV", "Wilson&Yang 2002"),
                      (5.0, 2e-12, "Chandra 5 keV", "Wilson&Yang 2002")],
        "Cyg X-1":   [(2.0, 3e-9,  "RXTE 2 keV",     "Wilms+2006"),
                      (10.0, 2e-9, "RXTE 10 keV",    "Wilms+2006"),
                      (50.0, 5e-10, "INTEGRAL 50 keV", "Wilms+2006")],
    }

    import astropy.units as u
    sed_obj = sedmod.SED(target_name=target_choice)
    for label, wl_um, fnu_jy, band, src in SEDS[target_choice]:
        sed_obj.add(sedmod.SEDPoint(
            label=label, wavelength=wl_um * u.micron,
            flux_density=fnu_jy * u.Jy, band=band, source=src,
        ))
    # X-ray points: provided as νF_ν in erg/s/cm² directly
    for E_keV, nfn, label, src in XRAY[target_choice]:
        sed_obj.add(sedmod.SEDPoint(
            label=label, energy=E_keV * u.keV,
            nu_f_nu=nfn * (u.erg / u.s / u.cm**2),
            band="xray", source=src,
        ))

    overlay = st.checkbox(
        "Overlay Elvis+1994 mean quasar template (schematic)",
        value=False,
        help="Hand-tabulated normalization — for visual context only."
    )

    fig = sedmod.render_sed(sed_obj, overplot_quasar_template=overlay)
    st.pyplot(fig, use_container_width=True)

    with st.expander("Data points and sources"):
        import pandas as pd
        rows = []
        for p in sed_obj.points:
            rows.append({
                "Label": p.label, "Band": p.band, "Source": p.source,
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
        st.caption(
            "Values are literature photometry from the sources noted, not "
            "extracted from the local FITS cutouts. Aperture photometry from "
            "the cutouts is planned for Phase 2."
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

            from blackhole.physics.variability import (fractional_rms,
                                                       fractional_rms_error)
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
