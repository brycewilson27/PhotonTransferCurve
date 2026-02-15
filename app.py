"""
Photon Transfer Curve (PTC) Analysis -- Streamlit Dashboard

Interactive sensor characterization tool for AstroTracker CMOS star tracker.
All analysis in native DN (stored 8-bit = native 12-bit, clamped at 255).
NO rescaling applied -- 1 stored DN = 1 native DN.

Supports two modes:
  - LIVE MODE: FITS data present locally, full interactive analysis
  - DEMO MODE: precomputed results from demo_data.json (Streamlit Cloud)
"""

from __future__ import annotations

import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from ptc_analysis import (
    ADU_MAX,
    DEFAULT_ROI,
    DEFAULT_SAT_THRESH,
    EMVA_SPEC,
    apply_quantization_correction,
)

# ---------------------------------------------------------------------------
# App configuration
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="PTC Sensor Analysis",
    page_icon=":telescope:",
    layout="wide",
)

BASE_DIR = Path(__file__).parent
ARTIFACTS = BASE_DIR / "PhotonTransferCurveArtifacts"
DEMO_DATA_PATH = BASE_DIR / "demo_data.json"

# Detect mode: LIVE if FITS data exists, DEMO otherwise
SWEEP_DIRS = {
    "Exposure Sweep (9 dB gain)": ARTIFACTS / "ExposureSweep_Gain9db",
    "Gain Sweep (10 ms exposure)": ARTIFACTS / "GainSweep_10msExpoTime",
}
GLOBAL_BIAS_PATHS = [
    ARTIFACTS / "ExposureSweep_Gain9db" / "BlackImage.fits",
    ARTIFACTS / "ExposureSweep_Gain9db" / "Black_Image2.fits",
]

_exposure_dir = ARTIFACTS / "ExposureSweep_Gain9db"
LIVE_MODE = _exposure_dir.is_dir() and any(_exposure_dir.glob("*/*.fits"))

# Image dimensions (from sensor spec)
IMG_HEIGHT = 1552
IMG_WIDTH = 2064


# ---------------------------------------------------------------------------
# Demo data loader (cached)
# ---------------------------------------------------------------------------
@st.cache_data
def load_demo_data() -> dict:
    """Load precomputed analysis results from demo_data.json."""
    with open(DEMO_DATA_PATH) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Conditional imports for live mode only
# ---------------------------------------------------------------------------
if LIVE_MODE:
    from astropy.io import fits

    from ptc_analysis import (
        analyze_gain_sweep_per_point,
        collect_flat_pairs,
        fit_gain_model,
        load_fits_image,
        make_master_bias,
        mean_var_from_pair,
        run_ptc_analysis,
    )


# ---------------------------------------------------------------------------
# Sidebar controls
# ---------------------------------------------------------------------------
def build_sidebar() -> dict:
    """Render sidebar controls and return current settings."""
    st.sidebar.title("PTC Analysis Controls")

    if not LIVE_MODE:
        st.sidebar.info(
            "Running in **demo mode** with precomputed results. "
            "Sidebar controls show default settings used during analysis."
        )

    dataset_name = st.sidebar.selectbox("Dataset", list(SWEEP_DIRS.keys()))
    is_gain_sweep = "Gain" in dataset_name

    st.sidebar.markdown("---")
    st.sidebar.subheader("Region of Interest")
    col1, col2 = st.sidebar.columns(2)
    roi_x = col1.number_input("Center X", 0, IMG_WIDTH - 1, DEFAULT_ROI[0])
    roi_y = col2.number_input("Center Y", 0, IMG_HEIGHT - 1, DEFAULT_ROI[1])
    col3, col4 = st.sidebar.columns(2)
    roi_w = col3.number_input("Width", 10, IMG_WIDTH, DEFAULT_ROI[2], step=10)
    roi_h = col4.number_input("Height", 10, IMG_HEIGHT, DEFAULT_ROI[3], step=10)
    roi = (int(roi_x), int(roi_y), int(roi_w), int(roi_h))

    st.sidebar.markdown("---")
    st.sidebar.subheader("Analysis Settings")
    sat_thresh = st.sidebar.slider(
        "Saturation threshold (%)",
        min_value=0.1,
        max_value=5.0,
        value=DEFAULT_SAT_THRESH * 100,
        step=0.1,
        help="Exclude pairs where more than this % of ROI pixels hit 0 or 255",
    ) / 100.0

    apply_quant = st.sidebar.checkbox("Quantization correction", value=True,
                                      help="Subtract 1/12 DN^2 quantization floor")

    st.sidebar.markdown("---")
    st.sidebar.subheader("Bias Method")
    if is_gain_sweep:
        st.sidebar.info("Gain Sweep has no global bias. Using local dark pairs.")
        use_local_dark = True
    else:
        use_local_dark = st.sidebar.radio(
            "Bias subtraction",
            ["Local dark pairs (per subfolder)", "Global master bias (BlackImage)"],
            index=0,
        ) == "Local dark pairs (per subfolder)"

    return {
        "dataset_name": dataset_name,
        "sweep_folder": SWEEP_DIRS[dataset_name],
        "is_gain_sweep": is_gain_sweep,
        "roi": roi,
        "sat_thresh": sat_thresh,
        "apply_quant": apply_quant,
        "use_local_dark": use_local_dark,
    }


# ---------------------------------------------------------------------------
# Analysis runner (live or demo)
# ---------------------------------------------------------------------------
@st.cache_data(show_spinner="Running PTC analysis...")
def cached_ptc_analysis(
    sweep_folder: str,
    bias_paths_str: list[str] | None,
    roi: tuple[int, int, int, int],
    sat_thresh: float,
    apply_quant: bool,
    use_local_dark: bool,
) -> tuple[dict, list[dict]]:
    """Run PTC analysis and return serializable dicts (for Streamlit caching)."""
    bias_paths = [Path(p) for p in bias_paths_str] if bias_paths_str else None
    fit_result, pair_results = run_ptc_analysis(
        sweep_folder=Path(sweep_folder),
        bias_paths=bias_paths,
        roi=roi,
        sat_thresh=sat_thresh,
        apply_quant_corr=apply_quant,
        use_local_dark=use_local_dark,
    )
    fit_dict = {
        "K_slope": fit_result.K_slope,
        "intercept": fit_result.intercept,
        "r_squared": fit_result.r_squared,
        "conversion_gain": fit_result.conversion_gain,
        "read_noise_dn": fit_result.read_noise_dn,
        "read_noise_e": fit_result.read_noise_e,
        "fwc_e": fit_result.fwc_e,
        "dynamic_range_db": fit_result.dynamic_range_db,
        "snr_max": fit_result.snr_max,
        "n_points_total": fit_result.n_points_total,
        "n_points_used": fit_result.n_points_used,
        "fit_mask": fit_result.fit_mask.tolist(),
    }
    pair_dicts = [
        {
            "label": pr.label,
            "mean_signal": pr.mean_signal,
            "variance": pr.variance,
            "sat_hi": pr.sat_hi,
            "sat_lo": pr.sat_lo,
            "path_a": str(pr.path_a),
            "path_b": str(pr.path_b),
        }
        for pr in pair_results
    ]
    return fit_dict, pair_dicts


def run_analysis(settings: dict) -> tuple[dict, list[dict]]:
    """Get analysis results -- live or from demo data."""
    if LIVE_MODE:
        bias_str = (
            [str(p) for p in GLOBAL_BIAS_PATHS]
            if GLOBAL_BIAS_PATHS[0].exists()
            else None
        )
        return cached_ptc_analysis(
            sweep_folder=str(settings["sweep_folder"]),
            bias_paths_str=bias_str,
            roi=settings["roi"],
            sat_thresh=settings["sat_thresh"],
            apply_quant=settings["apply_quant"],
            use_local_dark=settings["use_local_dark"],
        )
    else:
        demo = load_demo_data()
        key = "gain_sweep" if settings["is_gain_sweep"] else "exposure_sweep"
        return demo[key]["fit"], demo[key]["pairs"]


def get_demo_sweep(settings: dict) -> dict:
    """Get the full demo data dict for the selected sweep."""
    demo = load_demo_data()
    key = "gain_sweep" if settings["is_gain_sweep"] else "exposure_sweep"
    return demo[key]


@st.cache_data(show_spinner="Running per-point gain analysis...")
def cached_gain_per_point(
    sweep_folder: str,
    roi: tuple[int, int, int, int],
    sigma_read_e_anchor: float,
    apply_quant: bool,
    sat_thresh: float,
) -> tuple[list[dict], dict]:
    """Run per-point gain analysis and gain model fit, return serializable dicts."""
    gp_results = analyze_gain_sweep_per_point(
        sweep_folder=Path(sweep_folder),
        roi=roi,
        sigma_read_e_anchor=sigma_read_e_anchor,
        apply_quant_corr=apply_quant,
    )
    per_point = []
    for gp in gp_results:
        per_point.append({
            "gain_db": gp.gain_db,
            "mean_signal": float(gp.mean_signal),
            "pair_variance": float(gp.pair_variance),
            "K_apparent": float(gp.K_apparent) if np.isfinite(gp.K_apparent) else None,
            "K_true": float(gp.K_true) if np.isfinite(gp.K_true) else None,
            "dark_pair_variance": float(gp.dark_pair_variance),
            "read_noise_dn": float(gp.read_noise_dn),
            "read_noise_e": float(gp.read_noise_e) if np.isfinite(gp.read_noise_e) else None,
            "n_electrons": float(gp.n_electrons) if np.isfinite(gp.n_electrons) else None,
            "sat_hi": float(gp.sat_hi),
            "sat_lo": float(gp.sat_lo),
            "gain_mode": gp.gain_mode,
        })
    gm = fit_gain_model(gp_results, sat_thresh=sat_thresh)
    model_dict = {
        "K_0": float(gm.K_0),
        "dB_per_decade": float(gm.dB_per_decade),
        "dB_per_doubling": float(gm.dB_per_doubling),
        "scale_factor": float(gm.scale_factor),
        "r_squared": float(gm.r_squared),
        "n_points_used": gm.n_points_used,
        "gain_range": list(gm.gain_range),
        "sigma_read_e_mean": float(gm.sigma_read_e_mean),
    }
    return per_point, model_dict


def get_gain_analysis(settings: dict, exposure_fit: dict | None = None) -> tuple[list[dict], dict]:
    """Get per-point gain analysis and model -- live or demo."""
    if LIVE_MODE:
        # Use exposure sweep read_noise_e as anchor; fall back to 2.78
        anchor = 2.78
        if exposure_fit is not None:
            anchor = exposure_fit.get("read_noise_e", 2.78)
        return cached_gain_per_point(
            sweep_folder=str(settings["sweep_folder"]),
            roi=settings["roi"],
            sigma_read_e_anchor=anchor,
            apply_quant=settings["apply_quant"],
            sat_thresh=settings["sat_thresh"],
        )
    else:
        sweep_data = get_demo_sweep(settings)
        return sweep_data["per_point_analysis"], sweep_data["gain_model"]


# ---------------------------------------------------------------------------
# Helper: build a DataFrame from pair results
# ---------------------------------------------------------------------------
def pairs_to_dataframe(pair_dicts: list[dict], fit_dict: dict, apply_quant: bool) -> pd.DataFrame:
    """Build a summary DataFrame from pair result dicts."""
    mask = fit_dict["fit_mask"]
    variances_raw = [p["variance"] for p in pair_dicts]
    if apply_quant:
        variances = apply_quantization_correction(np.array(variances_raw)).tolist()
    else:
        variances = variances_raw

    rows = []
    for i, p in enumerate(pair_dicts):
        rows.append({
            "Label": p["label"],
            "Mean Signal (DN)": round(p["mean_signal"], 2),
            "Variance (DN^2)": round(variances[i], 2),
            "Sat Hi (%)": round(p["sat_hi"] * 100, 2),
            "Sat Lo (%)": round(p["sat_lo"] * 100, 2),
            "Used in Fit": mask[i],
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Page: Sensor Overview
# ---------------------------------------------------------------------------
def page_sensor_overview(settings: dict):
    st.header("Sensor Overview")

    if LIVE_MODE:
        sweep = settings["sweep_folder"]
        flat_pairs = collect_flat_pairs(sweep)
        n_pairs = len(flat_pairs)
        labels = [fp.label for fp in flat_pairs]
    else:
        sweep_data = get_demo_sweep(settings)
        n_pairs = sweep_data["n_pairs"]
        labels = sweep_data["labels"]

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Dataset Info")
        st.markdown(f"**Dataset:** {settings['dataset_name']}")
        if LIVE_MODE:
            st.markdown(f"**Path:** `{settings['sweep_folder']}`")
        st.markdown(f"**Exposure/gain levels:** {n_pairs}")
        st.markdown(f"**Frames per level:** 2 bright + 2 dark")
        st.markdown(f"**Total usable pairs:** {n_pairs}")

        if settings["is_gain_sweep"]:
            st.markdown(f"**Gain levels:** {', '.join(labels)}")
        else:
            st.markdown(f"**Exposure times:** {', '.join(labels)}")

    with col2:
        st.subheader("Sensor Specifications")
        st.markdown("**Sensor:** AstroTracker CMOS star tracker")
        st.markdown(f"**Underlying sensor:** {EMVA_SPEC['sensor_model']} ({EMVA_SPEC['sensor_type']})")
        st.markdown(f"**Pixel size:** {EMVA_SPEC['pixel_size_um']} um | **Diagonal:** {EMVA_SPEC['diagonal_mm']} mm")
        st.markdown("**Native ADC:** 12-bit (0-4095 DN)")
        st.markdown("**Stored format:** 8-bit via CLAMP (0-255 DN)")
        st.markdown(f"**Image size:** {EMVA_SPEC['resolution']} pixels")
        st.markdown("**Key:** 1 stored DN = 1 native DN (clamped, not rescaled)")

        if LIVE_MODE and flat_pairs:
            first_file = flat_pairs[0].bright_a
            try:
                with fits.open(str(first_file)) as hdul:
                    header = hdul[0].header
                    if len(header) > 0:
                        st.markdown("**FITS Header (sample):**")
                        header_items = []
                        for key in header.keys():
                            if key and key not in ("", "COMMENT", "HISTORY"):
                                header_items.append(f"  {key} = {header[key]}")
                        if header_items:
                            st.code("\n".join(header_items[:15]), language="text")
            except Exception:
                pass
        elif not LIVE_MODE and not settings["is_gain_sweep"]:
            sweep_data = get_demo_sweep(settings)
            if "header_text" in sweep_data:
                st.markdown("**FITS Header (sample):**")
                st.code(sweep_data["header_text"], language="text")

    # --- EMVA 1288 Reference Specifications ---
    st.markdown("---")
    st.subheader("EMVA 1288 Reference Specifications")
    st.markdown(
        f"The **{EMVA_SPEC['sensor_model']}** has been independently characterized by "
        f"**{EMVA_SPEC['test_lab']}** under the **{EMVA_SPEC['test_standard']}** standard "
        f"at {EMVA_SPEC['test_illumination']}, {EMVA_SPEC['test_temperature_C']}C. "
        "Two operating points are reported: **Low Conversion Gain (LCG)** at 0 dB analog gain, "
        "and **High Conversion Gain (HCG)** at 3 dB analog gain. Both use 12-bit ADC output."
    )

    lcg = EMVA_SPEC["LCG"]
    hcg = EMVA_SPEC["HCG"]

    spec_df = pd.DataFrame({
        "Parameter": [
            "Analog Gain",
            "Conversion Gain Mode",
            "System Gain K (DN/e-)",
            "Conversion Gain 1/K (e-/DN)",
            "Quantum Efficiency (%)",
            "Read Noise (e- rms)",
            "Read Noise (DN rms)",
            "Full Well Capacity (e-)",
            "SNR_max",
            "Dynamic Range (dB)",
            "Dynamic Range (bits)",
            "DSNU (e-)",
            "PRNU (%)",
            "Linearity Error (%)",
            "Black Level (DN)",
        ],
        "LCG (0 dB)": [
            "0 dB",
            "Low",
            f"{lcg['K_dn_per_e']:.3f}",
            f"{lcg['conv_gain_e_per_dn']:.3f}",
            f"{lcg['qe_percent']:.2f}",
            f"{lcg['read_noise_e']:.3f}",
            f"{lcg['read_noise_dn']:.3f}",
            f"{lcg['fwc_e']:,}",
            f"{lcg['snr_max']} ({lcg['snr_max_db']:.1f} dB)",
            f"{lcg['dynamic_range_db']:.1f}",
            f"{lcg['dynamic_range_bits']:.1f}",
            f"{lcg['dsnu_e']:.1f}",
            f"{lcg['prnu_percent']:.1f}",
            f"{lcg['linearity_error_min_pct']:.3f} to {lcg['linearity_error_max_pct']:.3f}",
            f"{lcg['black_level_dn']}",
        ],
        "HCG (3 dB)": [
            "3 dB",
            "High",
            f"{hcg['K_dn_per_e']:.3f}",
            f"{hcg['conv_gain_e_per_dn']:.3f}",
            f"{hcg['qe_percent']:.2f}",
            f"{hcg['read_noise_e']:.3f}",
            f"{hcg['read_noise_dn']:.3f}",
            f"{hcg['fwc_e']:,}",
            f"{hcg['snr_max']} ({hcg['snr_max_db']:.1f} dB)",
            f"{hcg['dynamic_range_db']:.1f}",
            f"{hcg['dynamic_range_bits']:.1f}",
            f"{hcg['dsnu_e']:.1f}",
            f"{hcg['prnu_percent']:.1f}",
            f"{hcg['linearity_error_min_pct']:.3f} to {hcg['linearity_error_max_pct']:.3f}",
            f"{hcg['black_level_dn']}",
        ],
    })
    st.dataframe(spec_df, use_container_width=True, hide_index=True)

    with st.expander("How do these specs relate to our measurements?"):
        st.markdown(
            "The EMVA 1288 specs characterize the **bare sensor** at low analog gain "
            "(0 or 3 dB) with **12-bit output**. Our AstroTracker operates at **9 dB "
            "analog gain** with **8-bit clamped output**, which creates key differences:\n\n"
            "- **System gain K** is much higher in our setup (3.13 DN/e- at 9 dB vs "
            "0.40 DN/e- at 0 dB) because analog gain amplifies the signal before digitization.\n"
            "- **FWC appears tiny** in our data (81 e-) because the 8-bit clamp at 255 DN "
            "cuts off the signal far before the sensor's true ~9,500 e- well fills. "
            "The true FWC is only accessible with 12-bit unclamped data.\n"
            "- **Read noise in electrons** is comparable (2.78 e- at 9 dB vs 1.95-5.56 e- "
            "in the spec), which makes physical sense -- read noise originates at the pixel "
            "level before amplification.\n"
            "- **Quantum efficiency** (~90%) is an intrinsic pixel property and does not depend "
            "on gain or bit depth. We cannot measure QE from our PTC (it requires a calibrated "
            "light source).\n\n"
            "The gain architecture has two layers: a **hardware conversion gain mode** "
            "(LCG/HCG, which changes pixel capacitance) and an **analog gain amplifier** "
            "(the dB setting). Our 9 dB setting applies amplification on top of "
            "whichever CG mode the AstroTracker firmware selects."
        )

    # Sample image with ROI overlay
    st.subheader("Sample Image with ROI Overlay")

    if LIVE_MODE and flat_pairs:
        sample_path = flat_pairs[len(flat_pairs) // 2].bright_a
        img = load_fits_image(sample_path)
        img_title = sample_path.name
    elif not LIVE_MODE and not settings["is_gain_sweep"]:
        sweep_data = get_demo_sweep(settings)
        if "sample_image_ds" in sweep_data:
            img = np.array(sweep_data["sample_image_ds"])
            ds = sweep_data["sample_image_ds_factor"]
            img_title = f"{sweep_data['sample_image_name']} (downsampled {ds}x)"
        else:
            img = None
            img_title = ""
    else:
        img = None
        img_title = ""

    if img is not None:
        roi = settings["roi"]
        rx, ry, rw, rh = roi
        # Scale ROI coords for downsampled image in demo mode
        if not LIVE_MODE and not settings["is_gain_sweep"]:
            ds = get_demo_sweep(settings).get("sample_image_ds_factor", 1)
            rx_s, ry_s, rw_s, rh_s = rx // ds, ry // ds, rw // ds, rh // ds
        else:
            rx_s, ry_s, rw_s, rh_s = rx, ry, rw, rh

        fig = go.Figure()
        fig.add_trace(go.Heatmap(
            z=img, colorscale="gray", showscale=True,
            colorbar=dict(title="DN"), reversescale=False,
        ))
        fig.add_shape(
            type="rect",
            x0=rx_s, y0=ry_s, x1=rx_s + rw_s, y1=ry_s + rh_s,
            line=dict(color="red", width=2),
        )
        fig.add_annotation(
            x=rx_s + rw_s / 2, y=ry_s - 10,
            text=f"ROI ({rw}x{rh})",
            showarrow=False,
            font=dict(color="red", size=12),
        )
        fig.update_layout(
            title=f"Sample: {img_title}",
            xaxis_title="X (pixels)", yaxis_title="Y (pixels)",
            yaxis=dict(autorange="reversed"),
            width=900, height=700,
        )
        st.plotly_chart(fig, use_container_width=True)
    elif not LIVE_MODE and settings["is_gain_sweep"]:
        st.info("Sample image preview is available for the Exposure Sweep dataset.")
    else:
        st.warning("No flat-field pairs found in this dataset.")


# ---------------------------------------------------------------------------
# Page: PTC Analysis
# ---------------------------------------------------------------------------
def page_ptc_analysis(settings: dict):
    st.header("Photon Transfer Curve Analysis")

    try:
        fit_dict, pair_dicts = run_analysis(settings)
    except ValueError as e:
        st.error(f"Analysis failed: {e}")
        return

    mask = np.array(fit_dict["fit_mask"])
    means = np.array([p["mean_signal"] for p in pair_dicts])
    variances_raw = np.array([p["variance"] for p in pair_dicts])
    if settings["apply_quant"]:
        variances = apply_quantization_correction(variances_raw)
    else:
        variances = variances_raw

    K = fit_dict["K_slope"]
    intercept = fit_dict["intercept"]

    # Fit line x range
    fit_x = np.linspace(0, min(means[mask].max() * 1.2, ADU_MAX), 200)
    fit_y = K * fit_x + intercept

    # --- Linear PTC plot ---
    use_log = st.checkbox("Log-log scale", value=False)

    fig = go.Figure()

    # Excluded points
    if (~mask).any():
        fig.add_trace(go.Scatter(
            x=means[~mask],
            y=variances[~mask],
            mode="markers",
            name="Excluded (saturated)",
            marker=dict(color="red", size=8, symbol="x"),
            text=[pair_dicts[i]["label"] for i in range(len(pair_dicts)) if not mask[i]],
            hovertemplate="<b>%{text}</b><br>Mean: %{x:.1f} DN<br>Var: %{y:.2f} DN^2<extra></extra>",
        ))

    # Used points
    fig.add_trace(go.Scatter(
        x=means[mask],
        y=variances[mask],
        mode="markers",
        name="Used in fit",
        marker=dict(color="dodgerblue", size=8),
        text=[pair_dicts[i]["label"] for i in range(len(pair_dicts)) if mask[i]],
        hovertemplate="<b>%{text}</b><br>Mean: %{x:.1f} DN<br>Var: %{y:.2f} DN^2<extra></extra>",
    ))

    # Fit line
    fig.add_trace(go.Scatter(
        x=fit_x, y=fit_y,
        mode="lines",
        name=f"Fit: Var = {K:.4f} * Mean + {intercept:.2f}",
        line=dict(color="orange", dash="dash"),
    ))

    # Vertical clamp ceiling line
    fig.add_vline(
        x=ADU_MAX, line_dash="dot", line_color="gray",
        annotation_text="255 DN clamp", annotation_position="top left",
    )

    axis_type = "log" if use_log else "linear"
    fig.update_layout(
        title="Photon Transfer Curve (Variance vs Mean Signal)",
        xaxis_title="Mean Signal (DN)",
        yaxis_title="Variance (DN^2)",
        xaxis=dict(type=axis_type),
        yaxis=dict(type=axis_type),
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
        height=550,
    )
    st.plotly_chart(fig, use_container_width=True)

    # Fit summary
    col1, col2, col3 = st.columns(3)
    col1.metric("K (DN/e-)", f"{K:.4f}")
    col2.metric("R-squared", f"{fit_dict['r_squared']:.4f}")
    col3.metric("Points used", f"{fit_dict['n_points_used']} / {fit_dict['n_points_total']}")

    # Derived sensor metrics
    st.markdown("---")
    st.subheader("Derived Sensor Metrics")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Conversion Gain", f"{fit_dict['conversion_gain']:.2f} e-/DN")
    m2.metric("Read Noise", f"{fit_dict['read_noise_e']:.2f} e- rms")
    m3.metric("Read Noise (DN)", f"{fit_dict['read_noise_dn']:.2f} DN rms")
    m4.metric("Dynamic Range", f"{fit_dict['dynamic_range_db']:.1f} dB")

    m5, m6, m7 = st.columns(3)
    m5.metric("Full Well Capacity", f"{fit_dict['fwc_e']:.0f} e-")
    m6.metric("SNR_max", f"{fit_dict['snr_max']:.1f}")
    m7.metric("Intercept (sigma^2_read)", f"{fit_dict['intercept']:.2f} DN^2")

    st.warning(
        "FWC is clamp-limited at 255 DN -- this is the capacity through the 8-bit output "
        "window, NOT the true sensor full well capacity. True FWC requires unclamped 12-bit "
        "data to measure."
    )

    # --- EMVA 1288 Spec Comparison ---
    if not settings["is_gain_sweep"]:
        st.markdown("---")
        st.subheader("Comparison with EMVA 1288 Spec Values")
        st.markdown(
            "Side-by-side comparison of our PTC-measured values (9 dB gain, 8-bit clamp) "
            "against the manufacturer EMVA 1288 specs (0-3 dB gain, 12-bit output)."
        )

        lcg = EMVA_SPEC["LCG"]
        hcg = EMVA_SPEC["HCG"]
        cmp_df = pd.DataFrame({
            "Metric": [
                "System Gain K (DN/e-)",
                "Conversion Gain (e-/DN)",
                "Read Noise (e- rms)",
                "Read Noise (DN rms)",
                "Full Well Capacity (e-)",
                "Dynamic Range (dB)",
                "SNR_max",
            ],
            "Our PTC (9 dB, 8-bit)": [
                f"{fit_dict['K_slope']:.4f}",
                f"{fit_dict['conversion_gain']:.3f}",
                f"{fit_dict['read_noise_e']:.2f}",
                f"{fit_dict['read_noise_dn']:.2f}",
                f"{fit_dict['fwc_e']:.0f} (clamp-limited)",
                f"{fit_dict['dynamic_range_db']:.1f}",
                f"{fit_dict['snr_max']:.1f}",
            ],
            "Spec LCG (0 dB, 12-bit)": [
                f"{lcg['K_dn_per_e']:.4f}",
                f"{lcg['conv_gain_e_per_dn']:.3f}",
                f"{lcg['read_noise_e']:.3f}",
                f"{lcg['read_noise_dn']:.3f}",
                f"{lcg['fwc_e']:,}",
                f"{lcg['dynamic_range_db']:.1f}",
                f"{lcg['snr_max']}",
            ],
            "Spec HCG (3 dB, 12-bit)": [
                f"{hcg['K_dn_per_e']:.4f}",
                f"{hcg['conv_gain_e_per_dn']:.3f}",
                f"{hcg['read_noise_e']:.3f}",
                f"{hcg['read_noise_dn']:.3f}",
                f"{hcg['fwc_e']:,}",
                f"{hcg['dynamic_range_db']:.1f}",
                f"{hcg['snr_max']}",
            ],
        })
        st.dataframe(cmp_df, use_container_width=True, hide_index=True)

        st.info(
            "**Why the big differences?** Our 9 dB analog gain amplifies the signal ~7.8x more than "
            "LCG mode (K=3.13 vs 0.40 DN/e-). The 8-bit clamp at 255 DN then limits the measurable "
            "FWC to just 81 e-. The spec's true FWC of ~9,500 e- (LCG) is only accessible "
            "with 12-bit unclamped output at low gain."
        )


# ---------------------------------------------------------------------------
# Page: Sensor Metrics Dashboard
# ---------------------------------------------------------------------------
def page_metrics_dashboard(settings: dict):
    st.header("Sensor Metrics Dashboard")

    try:
        fit_dict, pair_dicts = run_analysis(settings)
    except ValueError as e:
        st.error(f"Analysis failed: {e}")
        return

    # Metric cards
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Conversion Gain", f"{fit_dict['conversion_gain']:.2f} e-/DN")
    c2.metric("Read Noise", f"{fit_dict['read_noise_e']:.2f} e- rms")
    c3.metric("Dynamic Range", f"{fit_dict['dynamic_range_db']:.1f} dB")
    c4.metric("SNR_max", f"{fit_dict['snr_max']:.1f}")

    # FWC with caveat
    st.metric("Full Well Capacity", f"{fit_dict['fwc_e']:.0f} e-")
    st.warning(
        "FWC is clamp-limited at 255 DN -- this is the capacity through the 8-bit output "
        "window, NOT the true sensor full well capacity. True FWC requires unclamped 12-bit "
        "data to measure."
    )

    st.markdown("---")

    # Additional metrics
    c5, c6, c7, c8 = st.columns(4)
    c5.metric("K (slope)", f"{fit_dict['K_slope']:.4f} DN/e-")
    c6.metric("Intercept", f"{fit_dict['intercept']:.2f} DN^2")
    c7.metric("Read Noise (DN)", f"{fit_dict['read_noise_dn']:.2f} DN rms")
    c8.metric("R-squared", f"{fit_dict['r_squared']:.6f}")

    # --- EMVA 1288 Spec Comparison on Metrics Dashboard ---
    if not settings["is_gain_sweep"]:
        st.markdown("---")
        st.subheader("EMVA 1288 Spec Comparison")

        lcg = EMVA_SPEC["LCG"]
        hcg = EMVA_SPEC["HCG"]

        # Show measured vs spec as metric cards with deltas
        st.markdown("**Our PTC vs manufacturer specs** (deltas show our value relative to spec):")
        sc1, sc2, sc3 = st.columns(3)
        sc1.metric(
            "K (DN/e-)",
            f"{fit_dict['K_slope']:.3f}",
            delta=f"{fit_dict['K_slope'] - hcg['K_dn_per_e']:+.3f} vs HCG spec",
            delta_color="off",
        )
        sc2.metric(
            "Read Noise (e-)",
            f"{fit_dict['read_noise_e']:.2f}",
            delta=f"{fit_dict['read_noise_e'] - hcg['read_noise_e']:+.2f} vs HCG spec",
            delta_color="off",
        )
        sc3.metric(
            "FWC (e-)",
            f"{fit_dict['fwc_e']:.0f}",
            delta=f"vs {lcg['fwc_e']:,} e- true spec (LCG)",
            delta_color="off",
        )
        st.caption(
            "Deltas are informational -- direct comparison is not meaningful because our "
            "setup uses 9 dB analog gain + 8-bit clamp vs the spec's 0-3 dB gain + 12-bit output. "
            "See the Sensor Overview page for the full EMVA 1288 spec table."
        )

    st.markdown("---")

    # Per-pair table
    st.subheader("Per-Pair Summary")
    df = pairs_to_dataframe(pair_dicts, fit_dict, settings["apply_quant"])
    st.dataframe(
        df.style.map(
            lambda v: "background-color: #d4edda" if v else "background-color: #f8d7da",
            subset=["Used in Fit"],
        ),
        use_container_width=True,
        hide_index=True,
    )

    # CSV export
    csv_buf = io.StringIO()
    df.to_csv(csv_buf, index=False)
    st.download_button(
        label="Export CSV",
        data=csv_buf.getvalue(),
        file_name="ptc_results.csv",
        mime="text/csv",
    )


# ---------------------------------------------------------------------------
# Page: Gain Sweep Characterization
# ---------------------------------------------------------------------------
def page_gain_sweep(settings: dict):
    st.header("Gain Sweep Characterization")
    st.markdown(
        "How does the sensor respond across its analog gain range? "
        "This page shows signal, dark level, noise, and saturation "
        "behavior at fixed 10 ms exposure as gain sweeps from 3 to 22 dB."
    )

    if not settings["is_gain_sweep"]:
        st.info(
            "Select the **Gain Sweep (10 ms exposure)** dataset in the sidebar "
            "to use this page."
        )
        return

    try:
        fit_dict, pair_dicts = run_analysis(settings)
    except ValueError as e:
        st.error(f"Analysis failed: {e}")
        return

    # Extract numeric gain values and dark frame stats
    gain_nums = []
    for label in [p["label"] for p in pair_dicts]:
        try:
            gain_nums.append(float(label.replace("db", "").replace("dB", "")))
        except ValueError:
            gain_nums.append(0.0)

    means = [p["mean_signal"] for p in pair_dicts]
    variances_raw = np.array([p["variance"] for p in pair_dicts])
    if settings["apply_quant"]:
        variances = apply_quantization_correction(variances_raw).tolist()
    else:
        variances = variances_raw.tolist()

    # Get dark frame stats from demo data or live analysis
    dark_means = []
    dark_stds = []
    if LIVE_MODE:
        flat_pairs = collect_flat_pairs(settings["sweep_folder"])
        roi = settings["roi"]
        rx, ry, rw, rh = roi
        for fp in flat_pairs:
            if fp.dark_a is not None:
                dark_img = load_fits_image(fp.dark_a)
                dark_roi = dark_img[ry:ry + rh, rx:rx + rw]
                dark_means.append(float(dark_roi.mean()))
                dark_stds.append(float(dark_roi.std()))
            else:
                dark_means.append(np.nan)
                dark_stds.append(np.nan)
    else:
        sweep_data = get_demo_sweep(settings)
        fe_data = sweep_data.get("frame_explorer", [])
        for fe in fe_data:
            ds = fe.get("dark_stats")
            if ds:
                dark_means.append(ds["ROI Mean"])
                dark_stds.append(ds["ROI Std"])
            else:
                dark_means.append(np.nan)
                dark_stds.append(np.nan)

    # --- Load per-point gain analysis ---
    # Get exposure sweep fit for the sigma_read_e anchor
    exposure_fit = None
    if LIVE_MODE:
        try:
            exp_settings = dict(settings, is_gain_sweep=False,
                                sweep_folder=SWEEP_DIRS["Exposure Sweep (9 dB gain)"])
            exposure_fit, _ = run_analysis(exp_settings)
        except Exception:
            pass
    else:
        try:
            demo = load_demo_data()
            exposure_fit = demo["exposure_sweep"]["fit"]
        except Exception:
            pass

    try:
        per_point, gain_model = get_gain_analysis(settings, exposure_fit)
    except Exception as exc:
        st.warning(f"Per-point gain analysis unavailable: {exc}")
        per_point, gain_model = [], {}

    # --- Overview metrics ---
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Gain Levels", f"{len(pair_dicts)}")
    c2.metric("Gain Range", f"{gain_nums[0]} -- {gain_nums[-1]} dB")
    c3.metric("Fixed Exposure", "10 ms")
    unsaturated = sum(1 for p in pair_dicts if p["sat_hi"] < settings["sat_thresh"])
    c4.metric("Unsaturated Levels", f"{unsaturated} / {len(pair_dicts)}")

    # --- 1. Signal Amplification ---
    st.markdown("---")
    st.subheader("Signal Amplification")
    st.markdown(
        "Higher analog gain amplifies the signal from a fixed photon count. "
        "At high gains the signal clips at 255 DN (8-bit clamp ceiling)."
    )

    fig_signal = go.Figure()
    fig_signal.add_trace(go.Scatter(
        x=gain_nums, y=means,
        mode="markers+lines", name="Bright (bias-subtracted)",
        marker=dict(size=7, color="dodgerblue"),
        text=[p["label"] for p in pair_dicts],
        hovertemplate="<b>%{text}</b><br>Mean: %{y:.1f} DN<extra></extra>",
    ))
    if dark_means:
        fig_signal.add_trace(go.Scatter(
            x=gain_nums, y=dark_means,
            mode="markers+lines", name="Dark frame (ROI)",
            marker=dict(size=6, color="gray", symbol="diamond"),
            hovertemplate="Gain: %{x} dB<br>Dark mean: %{y:.1f} DN<extra></extra>",
        ))
    fig_signal.add_hline(y=ADU_MAX, line_dash="dot", line_color="red",
                         annotation_text="255 DN clamp")
    fig_signal.update_layout(
        title="Mean Signal vs Analog Gain",
        xaxis_title="Analog Gain (dB)", yaxis_title="Mean Signal (DN)",
        height=450,
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
    )
    st.plotly_chart(fig_signal, use_container_width=True)

    # === NEW GAIN ANALYSIS SECTIONS ===
    if per_point and gain_model:
        sat_thresh_val = settings["sat_thresh"]

        # Filter to unsaturated points for plotting
        pp_unsat = [p for p in per_point if p["sat_hi"] < sat_thresh_val
                    and p.get("K_true") is not None]
        pp_all_valid = [p for p in per_point if p.get("K_true") is not None]

        # --- Section: System Gain (K) vs Analog Gain ---
        st.markdown("---")
        st.subheader("System Gain (K) vs Analog Gain")
        st.markdown(
            "Per-gain-level system gain estimates. **K_true** corrects for read-noise bias "
            "using the quadratic PTC equation; **K_apparent** (naive variance/mean) is biased high."
        )

        fig_k = go.Figure()

        # K_apparent (gray, all valid unsaturated)
        fig_k.add_trace(go.Scatter(
            x=[p["gain_db"] for p in pp_unsat],
            y=[p["K_apparent"] for p in pp_unsat],
            mode="markers",
            name="K_apparent (biased)",
            marker=dict(color="gray", size=7, symbol="diamond"),
            hovertemplate="Gain: %{x} dB<br>K_app: %{y:.3f} DN/e-<extra></extra>",
        ))

        # K_true (blue, unsaturated)
        fig_k.add_trace(go.Scatter(
            x=[p["gain_db"] for p in pp_unsat],
            y=[p["K_true"] for p in pp_unsat],
            mode="markers",
            name="K_true (corrected)",
            marker=dict(color="dodgerblue", size=8),
            hovertemplate="Gain: %{x} dB<br>K_true: %{y:.3f} DN/e-<extra></extra>",
        ))

        # Exponential model fit line
        gm = gain_model
        g_range = np.linspace(gm["gain_range"][0] - 5, gm["gain_range"][1] + 5, 200)
        K_model = gm["K_0"] * 10.0 ** (g_range / gm["dB_per_decade"])
        fig_k.add_trace(go.Scatter(
            x=g_range.tolist(), y=K_model.tolist(),
            mode="lines",
            name=f"Model: K = {gm['K_0']:.2f} * 10^(g/{gm['dB_per_decade']:.0f})",
            line=dict(color="orange", dash="dash", width=2),
        ))

        # ExposureSweep K reference star
        exp_K = exposure_fit["K_slope"] if exposure_fit else 3.13
        fig_k.add_trace(go.Scatter(
            x=[9], y=[exp_K],
            mode="markers",
            name=f"ExposureSweep K = {exp_K:.2f}",
            marker=dict(color="gold", size=14, symbol="star",
                        line=dict(color="black", width=1.5)),
        ))

        fig_k.update_layout(
            title="System Gain K vs Analog Gain",
            xaxis_title="Analog Gain (dB)",
            yaxis_title="K (DN/e-)",
            yaxis=dict(type="log"),
            height=500,
            legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
        )
        st.plotly_chart(fig_k, use_container_width=True)

        # Metric cards
        mc1, mc2, mc3, mc4 = st.columns(4)
        mc1.metric("K_0 (at 0 dB)", f"{gm['K_0']:.4f} DN/e-")
        mc2.metric("dB per decade", f"{gm['dB_per_decade']:.1f}")
        mc3.metric("dB per doubling", f"{gm['dB_per_doubling']:.1f}")
        mc4.metric("Model R^2", f"{gm['r_squared']:.4f}")

        with st.expander("How is K_true computed?"):
            st.markdown(
                "**K_apparent** = variance / mean overestimates K because the variance "
                "includes read noise. At 9 dB, this overestimate is ~28%.\n\n"
                "**K_true** solves the full single-point PTC equation:\n"
            )
            st.latex(r"\mathrm{Var} = K \cdot \bar{S} + K^2 \cdot \sigma_{\mathrm{read},e^-}^2")
            st.markdown("Rearranging as a quadratic in K:")
            st.latex(r"K^2 \cdot \sigma_{\mathrm{read},e^-}^2 + K \cdot \bar{S} - \mathrm{Var} = 0")
            st.markdown("Taking the positive root:")
            st.latex(
                r"K_{\mathrm{true}} = \frac{-\bar{S} + \sqrt{\bar{S}^2 + "
                r"4\,\sigma_{\mathrm{read},e^-}^2 \cdot \mathrm{Var}}}"
                r"{2\,\sigma_{\mathrm{read},e^-}^2}"
            )
            anchor_val = exposure_fit.get("read_noise_e", 2.78) if exposure_fit else 2.78
            st.markdown(
                f"The anchor value **sigma_read_e = {anchor_val:.2f} e-** comes from the "
                "ExposureSweep PTC y-intercept (the read noise from the full multi-point fit)."
            )

        # --- Section: Read Noise vs Analog Gain ---
        st.markdown("---")
        st.subheader("Read Noise vs Analog Gain")
        st.markdown(
            "Read noise measured from dark pair differences. In DN, it grows with gain "
            "(amplified fixed noise). In electrons (input-referred), it stays flat -- "
            "showing the constant noise floor of the readout electronics."
        )

        fig_rn = go.Figure()

        # Read noise in DN (growing)
        rn_gains = [p["gain_db"] for p in pp_unsat]
        rn_dn = [p["read_noise_dn"] for p in pp_unsat]
        rn_e = [p["read_noise_e"] for p in pp_unsat]

        fig_rn.add_trace(go.Scatter(
            x=rn_gains, y=rn_dn,
            mode="markers+lines", name="Read Noise (DN)",
            marker=dict(color="steelblue", size=7),
            hovertemplate="Gain: %{x} dB<br>RN: %{y:.2f} DN<extra></extra>",
        ))
        fig_rn.add_trace(go.Scatter(
            x=rn_gains, y=rn_e,
            mode="markers+lines", name="Read Noise (e-)",
            marker=dict(color="coral", size=7, symbol="triangle-up"),
            line=dict(dash="dash"),
            hovertemplate="Gain: %{x} dB<br>RN: %{y:.2f} e-<extra></extra>",
        ))

        # Mean sigma_read_e annotation
        fig_rn.add_hline(
            y=gm["sigma_read_e_mean"], line_dash="dot", line_color="coral",
            annotation_text=f"mean = {gm['sigma_read_e_mean']:.2f} e-",
            annotation_position="top right",
        )

        fig_rn.update_layout(
            title="Read Noise vs Analog Gain (Dark Pair-Difference Method)",
            xaxis_title="Analog Gain (dB)",
            yaxis_title="Read Noise",
            height=450,
            legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
        )
        st.plotly_chart(fig_rn, use_container_width=True)

        with st.expander("Why is read noise in electrons constant?"):
            st.markdown(
                "Read noise originates in the readout electronics **before** the analog gain "
                "amplifier. When measured in electrons (input-referred), it reflects the "
                "intrinsic noise floor of the pixel and is independent of gain.\n\n"
                "In DN (output-referred), read noise grows with gain because the amplifier "
                "scales both signal and noise equally.\n\n"
                f"**Dark pair-diff read noise:** {gm['sigma_read_e_mean']:.2f} e- "
                f"(mean across {gm['n_points_used']} unsaturated gains)\n\n"
                f"**PTC intercept read noise:** "
                f"{exposure_fit.get('read_noise_e', 2.78) if exposure_fit else 2.78:.2f} e- "
                "(from ExposureSweep fit)\n\n"
                "The dark pair-diff value is lower because it measures only the temporal "
                "readout noise floor, while the PTC intercept includes additional constant "
                "noise sources (e.g., low-level PRNU residuals, ADC non-linearity)."
            )

        # --- Section: Electron Count Consistency Check ---
        st.markdown("---")
        st.subheader("Electron Count Consistency Check")
        st.markdown(
            "With fixed illumination (10 ms), every gain level should collect the same "
            "number of photoelectrons. This is the key self-consistency check."
        )

        ne_gains = [p["gain_db"] for p in pp_unsat if p.get("n_electrons") is not None]
        ne_vals = [p["n_electrons"] for p in pp_unsat if p.get("n_electrons") is not None]

        if ne_vals:
            ne_mean = float(np.mean(ne_vals))
            ne_std = float(np.std(ne_vals))
            ne_cv = 100.0 * ne_std / ne_mean if ne_mean > 0 else 0

            fig_ne = go.Figure()
            fig_ne.add_trace(go.Scatter(
                x=ne_gains, y=ne_vals,
                mode="markers",
                name="n_electrons",
                marker=dict(color="mediumseagreen", size=9),
                hovertemplate="Gain: %{x} dB<br>n_e: %{y:.1f} e-<extra></extra>",
            ))

            # Mean +/- std band
            fig_ne.add_hrect(
                y0=ne_mean - ne_std, y1=ne_mean + ne_std,
                fillcolor="mediumseagreen", opacity=0.15,
                line_width=0,
                annotation_text=f"mean = {ne_mean:.1f} +/- {ne_std:.1f} e-",
                annotation_position="top right",
            )
            fig_ne.add_hline(y=ne_mean, line_dash="dash", line_color="green")

            fig_ne.update_layout(
                title="Electron Count vs Analog Gain (Unsaturated Points Only)",
                xaxis_title="Analog Gain (dB)",
                yaxis_title="Electrons (e-)",
                height=420,
            )
            st.plotly_chart(fig_ne, use_container_width=True)

            nc1, nc2, nc3 = st.columns(3)
            nc1.metric("Mean Electron Count", f"{ne_mean:.1f} e-")
            nc2.metric("Std Dev", f"{ne_std:.1f} e-")
            nc3.metric("CV", f"{ne_cv:.1f}%")

            st.markdown(
                "A low CV (< 10%) confirms that the gain correction is working correctly "
                "and the illumination is indeed constant across the sweep. Saturated points "
                "(> 15 dB) are excluded because clamp compression inflates the apparent electron count."
            )

        # --- Section: dB Scale Calibration ---
        st.markdown("---")
        st.subheader("dB Scale Calibration")
        st.markdown(
            "This table maps gain dB settings "
            "to equivalent standard voltage dB and the predicted system gain K."
        )

        cal_rows = []
        for g_db in [0, 3, 5, 7, 9, 11, 13, 15]:
            std_db = g_db / gm["scale_factor"]
            K_pred = gm["K_0"] * 10.0 ** (g_db / gm["dB_per_decade"])
            cal_rows.append({
                "Gain (dB)": g_db,
                "Std Voltage (dB)": round(std_db, 1),
                "Predicted K (DN/e-)": round(K_pred, 4),
            })

        # Add EMVA spec reference rows
        cal_rows.append({
            "Gain (dB)": "EMVA LCG (0 dB)",
            "Std Voltage (dB)": 0.0,
            "Predicted K (DN/e-)": EMVA_SPEC["LCG"]["K_dn_per_e"],
        })
        cal_rows.append({
            "Gain (dB)": "EMVA HCG (3 dB)",
            "Std Voltage (dB)": 3.0,
            "Predicted K (DN/e-)": EMVA_SPEC["HCG"]["K_dn_per_e"],
        })

        cal_df = pd.DataFrame(cal_rows)
        st.dataframe(cal_df, use_container_width=True, hide_index=True)

        st.markdown(
            f"**K_0 = {gm['K_0']:.2f} DN/e-** falls between EMVA LCG "
            f"({EMVA_SPEC['LCG']['K_dn_per_e']:.3f}) and HCG "
            f"({EMVA_SPEC['HCG']['K_dn_per_e']:.3f}), with ratios "
            f"{gm['K_0'] / EMVA_SPEC['LCG']['K_dn_per_e']:.2f}x and "
            f"{gm['K_0'] / EMVA_SPEC['HCG']['K_dn_per_e']:.2f}x respectively."
        )

        with st.expander("Are these standard dB?"):
            st.markdown(
                f"The gain values are close to standard voltage dB, with a small "
                f"deviation:\n\n"
                f"- Standard voltage dB: **6.02 dB** doubles the voltage (and K)\n"
                f"- Measured: **{gm['dB_per_doubling']:.1f} dB** doubles K\n"
                f"- Scale factor: **{gm['scale_factor']:.2f}x** "
                f"(measured / standard voltage dB)\n\n"
                f"At 9 dB gain, the actual amplification is a factor of "
                f"~{10.0 ** (9 / gm['dB_per_decade']):.1f}x "
                f"(standard dB would predict ~{10.0 ** (9 / 20.0):.1f}x)."
            )

    # --- 2. Dark Level (Black Level) ---
    if dark_means and not all(np.isnan(v) for v in dark_means):
        st.markdown("---")
        st.subheader("Dark Level (Black Level) vs Gain")
        st.markdown(
            "The dark frame captures the sensor's output with no illumination at each gain. "
            "As gain increases, thermal noise and readout offsets are amplified, "
            "raising the black level floor."
        )

        fig_dark = go.Figure()
        fig_dark.add_trace(go.Scatter(
            x=gain_nums, y=dark_means,
            mode="markers+lines", name="Dark ROI Mean",
            marker=dict(size=7, color="slategray"),
            line=dict(color="slategray"),
        ))
        fig_dark.add_trace(go.Scatter(
            x=gain_nums, y=dark_stds,
            mode="markers+lines", name="Dark ROI Std Dev",
            marker=dict(size=6, color="orange", symbol="triangle-up"),
            line=dict(color="orange", dash="dash"),
        ))
        fig_dark.update_layout(
            title="Dark Frame Statistics vs Analog Gain",
            xaxis_title="Analog Gain (dB)", yaxis_title="DN",
            height=420,
            legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
        )
        st.plotly_chart(fig_dark, use_container_width=True)

    # --- 3. Temporal Noise ---
    st.markdown("---")
    st.subheader("Temporal Noise vs Gain")
    st.markdown(
        "Temporal noise (frame-to-frame variance from the pair-difference method) "
        "increases with gain as both shot noise and read noise are amplified. "
        "At very high gains, variance drops because most pixels are clamped at 255."
    )

    noise_stds = [v ** 0.5 for v in variances]

    fig_noise = go.Figure()
    fig_noise.add_trace(go.Scatter(
        x=gain_nums, y=noise_stds,
        mode="markers+lines", name="Temporal Noise (std dev)",
        marker=dict(size=7, color="mediumseagreen"),
        text=[p["label"] for p in pair_dicts],
        hovertemplate="<b>%{text}</b><br>Noise: %{y:.2f} DN rms<extra></extra>",
    ))
    fig_noise.update_layout(
        title="Temporal Noise (Pair-Difference Std Dev) vs Analog Gain",
        xaxis_title="Analog Gain (dB)", yaxis_title="Noise (DN rms)",
        height=420,
    )
    st.plotly_chart(fig_noise, use_container_width=True)

    with st.expander("Why does noise drop at high gain?"):
        st.markdown(
            "At high gain, pixels are driven to the 255 DN clamp ceiling. "
            "Clamped pixels have zero temporal variation (they're always 255), "
            "so the measured variance **decreases** even though the underlying "
            "analog noise is still increasing. This is an artifact of the 8-bit "
            "output window, not a real noise reduction."
        )

    # --- 4. Saturation ---
    st.markdown("---")
    st.subheader("Saturation vs Gain")
    st.markdown(
        "Fraction of ROI pixels hitting the floor (0 DN) or ceiling (255 DN). "
        "Above ~15 dB, significant clipping begins."
    )

    sat_his = [p["sat_hi"] * 100 for p in pair_dicts]
    sat_los = [p["sat_lo"] * 100 for p in pair_dicts]

    fig_sat = go.Figure()
    fig_sat.add_trace(go.Scatter(
        x=gain_nums, y=sat_his,
        mode="markers+lines", name="Ceiling (255 DN)",
        marker=dict(color="red", size=6),
    ))
    fig_sat.add_trace(go.Scatter(
        x=gain_nums, y=sat_los,
        mode="markers+lines", name="Floor (0 DN)",
        marker=dict(color="blue", size=6),
    ))
    fig_sat.add_hline(y=settings["sat_thresh"] * 100, line_dash="dash",
                      line_color="orange", annotation_text="Exclusion threshold")
    fig_sat.update_layout(
        title="Saturation Fraction vs Analog Gain",
        xaxis_title="Analog Gain (dB)", yaxis_title="Saturated Pixels (%)",
        height=420,
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
    )
    st.plotly_chart(fig_sat, use_container_width=True)

    # --- 5. Signal-to-Noise Ratio ---
    st.markdown("---")
    st.subheader("Signal-to-Noise Ratio vs Gain")
    st.markdown(
        "SNR = mean signal / temporal noise. Initially SNR is relatively flat "
        "(both signal and noise scale with gain), but it degrades as saturation "
        "compresses the signal while the noise structure changes."
    )

    snr_vals = []
    for m, ns in zip(means, noise_stds):
        if ns > 0 and m > 0:
            snr_vals.append(m / ns)
        else:
            snr_vals.append(0.0)

    fig_snr = go.Figure()
    fig_snr.add_trace(go.Scatter(
        x=gain_nums, y=snr_vals,
        mode="markers+lines", name="SNR",
        marker=dict(size=7, color="purple"),
        text=[p["label"] for p in pair_dicts],
        hovertemplate="<b>%{text}</b><br>SNR: %{y:.1f}<extra></extra>",
    ))
    fig_snr.update_layout(
        title="SNR (Mean / Temporal Noise) vs Analog Gain",
        xaxis_title="Analog Gain (dB)", yaxis_title="SNR",
        height=420,
    )
    st.plotly_chart(fig_snr, use_container_width=True)

    # --- 6. Summary Table ---
    st.markdown("---")
    st.subheader("Per-Gain Summary Table")

    # Build lookup from per-point analysis by gain_db
    pp_lookup = {}
    if per_point:
        for pp in per_point:
            pp_lookup[pp["gain_db"]] = pp

    rows = []
    for i, p in enumerate(pair_dicts):
        row = {
            "Gain (dB)": gain_nums[i],
            "Mean Signal (DN)": round(means[i], 2),
            "Temporal Noise (DN)": round(noise_stds[i], 2),
            "SNR": round(snr_vals[i], 1),
            "Sat Ceiling (%)": round(sat_his[i], 3),
            "Sat Floor (%)": round(sat_los[i], 3),
        }
        if dark_means:
            row["Dark Mean (DN)"] = round(dark_means[i], 2) if not np.isnan(dark_means[i]) else "N/A"
            row["Dark Std (DN)"] = round(dark_stds[i], 2) if not np.isnan(dark_stds[i]) else "N/A"
        # Per-point gain analysis columns
        pp = pp_lookup.get(gain_nums[i])
        if pp:
            row["K_apparent (DN/e-)"] = round(pp["K_apparent"], 3) if pp.get("K_apparent") is not None else "-"
            row["K_true (DN/e-)"] = round(pp["K_true"], 3) if pp.get("K_true") is not None else "-"
            row["Read Noise (DN)"] = round(pp["read_noise_dn"], 2)
            row["Read Noise (e-)"] = round(pp["read_noise_e"], 2) if pp.get("read_noise_e") is not None else "-"
            row["n_electrons"] = round(pp["n_electrons"], 1) if pp.get("n_electrons") is not None else "-"
        rows.append(row)

    df = pd.DataFrame(rows)
    st.dataframe(df, use_container_width=True, hide_index=True)

    # CSV export
    csv_buf = io.StringIO()
    df.to_csv(csv_buf, index=False)
    st.download_button(
        label="Export CSV",
        data=csv_buf.getvalue(),
        file_name="gain_sweep_characterization.csv",
        mime="text/csv",
    )


# ---------------------------------------------------------------------------
# Page: PTC Theory & Derivation
# ---------------------------------------------------------------------------
def page_ptc_derivation(settings: dict):
    st.header("PTC Theory & Derivation")
    st.markdown(
        "A complete mathematical walkthrough of the Photon Transfer Curve method, "
        "including derivations for both exposure-varying and gain-varying sweep approaches."
    )

    # ---- 1. Introduction ----
    st.subheader("1. What Is a Photon Transfer Curve?")
    st.markdown(
        "The **Photon Transfer Curve (PTC)** is a fundamental tool for characterizing "
        "image sensors. It plots **temporal noise variance** against **mean signal** "
        "across a range of illumination levels. Because photon arrival is governed by "
        "Poisson statistics, the shape of this curve directly reveals the sensor's "
        "gain, read noise, and full-well capacity."
    )
    st.info(
        "The PTC is essentially a sensor's *fingerprint* -- a single experiment that "
        "extracts multiple performance parameters from first principles."
    )

    # ---- 2. Signal Chain Model ----
    st.subheader("2. Signal Chain Model")
    st.markdown(
        "Light incident on a pixel generates photoelectrons. The electrons are converted "
        "to voltage at the pixel's sense node, amplified by the **analog gain** stage, "
        "then digitized by the ADC into Digital Numbers (DN). The system gain **K** "
        "captures the entire chain:"
    )
    st.latex(r"\mathrm{Photons} \xrightarrow{\mathrm{QE}} e^- \xrightarrow{\mathrm{pixel}} V "
             r"\xrightarrow{\mathrm{gain}} V' \xrightarrow{\mathrm{ADC}} \mathrm{DN}")
    st.latex(r"S_{\mathrm{DN}} = K \cdot S_{e^-}")
    st.markdown("where:")
    st.markdown(
        "- $S_{\\mathrm{DN}}$ = signal in DN (digital numbers)\n"
        "- $S_{e^-}$ = signal in electrons\n"
        "- $K$ = system gain (DN/e$^-$), the PTC slope -- **includes analog gain**\n"
        "- Conversion gain $g = 1/K$ (e$^-$/DN), the Janesick convention"
    )

    with st.expander("Gain architecture of the IMX900"):
        st.markdown(
            "The **Sony IMX900** sensor has a **two-layer gain architecture**:\n\n"
            "1. **Conversion Gain Mode (LCG/HCG)** -- A hardware switch at the pixel "
            "that changes the capacitance of the sense node (charge-to-voltage conversion). "
            "HCG mode produces ~4.2x more DN per electron than LCG at the same analog gain.\n\n"
            "2. **Analog Gain (dB)** -- A programmable amplifier that scales the voltage "
            "before the ADC. Higher dB = more amplification = more DN per electron.\n\n"
            "The PTC-measured **K** bundles both layers together. At the spec's 0 dB (LCG), "
            f"K = {EMVA_SPEC['LCG']['K_dn_per_e']} DN/e-. At 3 dB (HCG), "
            f"K = {EMVA_SPEC['HCG']['K_dn_per_e']} DN/e-. Our AstroTracker "
            "at 9 dB measures K = 3.13 DN/e-, reflecting the amplification.\n\n"
            "**For PTC analysis and simulation, K is all you need** -- you don't need to "
            "know what fraction comes from conversion gain vs analog gain."
        )

    st.markdown("**Noise sources in the signal chain:**")
    st.markdown(
        "| Source | Domain | Variance |\n"
        "|--------|--------|----------|\n"
        "| Photon shot noise | electrons | $\\sigma^2_{\\mathrm{shot}} = S_{e^-}$ (Poisson) |\n"
        "| Read noise | electrons | $\\sigma^2_{\\mathrm{read}}$ (constant) |\n"
        "| Quantization noise | DN | $\\sigma^2_{q} = \\Delta^2 / 12$ |\n"
        "| Fixed-pattern noise | DN | $\\sigma^2_{\\mathrm{FPN}}$ (spatial, removed by pair method) |"
    )

    # ---- 3. Core PTC Derivation ----
    st.subheader("3. Core PTC Derivation")

    st.markdown("**Step 1: Photon statistics (Poisson process)**")
    st.markdown(
        "Photon arrival follows a Poisson distribution. For a Poisson random variable, "
        "the variance equals the mean:"
    )
    st.latex(r"\mathrm{Var}(S_{e^-}) = S_{e^-}")

    st.markdown("**Step 2: Transform to DN domain**")
    st.markdown("Since $S_{\\mathrm{DN}} = K \\cdot S_{e^-}$, we can write $S_{e^-} = S_{\\mathrm{DN}} / K$ and apply the variance scaling rule:")
    st.latex(r"\mathrm{Var}(S_{\mathrm{DN}}) = K^2 \cdot \mathrm{Var}(S_{e^-}) = K^2 \cdot S_{e^-} = K^2 \cdot \frac{S_{\mathrm{DN}}}{K} = K \cdot S_{\mathrm{DN}}")

    st.markdown("**Step 3: Add read noise**")
    st.markdown("Read noise is independent of signal and adds in quadrature (variances add):")
    st.latex(r"\boxed{\mathrm{Var}(S_{\mathrm{DN}}) = K \cdot S_{\mathrm{DN}} + \sigma^2_{\mathrm{read,DN}}}")
    st.success(
        "This is THE photon transfer equation. It is linear in mean signal: "
        "slope = K (DN/e-), y-intercept = read noise variance (DN^2)."
    )

    st.markdown("**Step 4: Quantization correction (optional)**")
    st.markdown(
        "The ADC introduces a uniform quantization noise with variance $\\Delta^2/12$, "
        "where $\\Delta$ is the step size. For our 12-bit ADC with 1 DN steps:"
    )
    st.latex(r"\sigma^2_q = \frac{1}{12} \approx 0.083 \;\mathrm{DN}^2")
    st.markdown("The corrected variance is:")
    st.latex(r"\mathrm{Var}_{\mathrm{corrected}} = \mathrm{Var}_{\mathrm{measured}} - \frac{1}{12}")

    # ---- 4. Pair-Difference Method ----
    st.subheader("4. Pair-Difference Method")
    st.markdown(
        "Fixed-pattern noise (FPN) is **spatial** variation that is constant between frames. "
        "It includes PRNU (photo-response non-uniformity) and DSNU (dark signal non-uniformity). "
        "If we compute variance from a single frame, FPN inflates the result. "
        "The pair-difference method eliminates it."
    )

    st.markdown("**Derivation:**")
    st.markdown(
        "Given two frames A and B taken under identical conditions, FPN is identical in both:"
    )
    st.latex(r"A = S + N_{\mathrm{temporal},A} + N_{\mathrm{FPN}}")
    st.latex(r"B = S + N_{\mathrm{temporal},B} + N_{\mathrm{FPN}}")
    st.markdown("Taking the difference cancels the FPN and the signal:")
    st.latex(r"A - B = N_{\mathrm{temporal},A} - N_{\mathrm{temporal},B}")
    st.markdown(
        "Since the temporal noise in A and B is independent and identically distributed:"
    )
    st.latex(r"\mathrm{Var}(A - B) = \mathrm{Var}(N_A) + \mathrm{Var}(N_B) = 2 \cdot \mathrm{Var}(N_{\mathrm{temporal}})")
    st.markdown("Therefore:")
    st.latex(r"\boxed{\mathrm{Var}_{\mathrm{temporal}} = \frac{1}{2} \mathrm{Var}(A - B)}")

    st.markdown("**Mean signal** is computed from the average of both frames after bias subtraction:")
    st.latex(r"\bar{S} = \frac{1}{2}(\bar{A} + \bar{B}) - \bar{B}_{\mathrm{bias}}")

    # ---- 5. Derived Metrics ----
    st.subheader("5. Derived Sensor Metrics")

    st.markdown("From the linear fit $\\mathrm{Var} = K \\cdot \\bar{S} + b$, we extract:")

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**System Gain**")
        st.latex(r"K = \text{slope} \quad [\mathrm{DN}/e^-]")

        st.markdown("**Conversion Gain**")
        st.latex(r"g = \frac{1}{K} \quad [e^-/\mathrm{DN}]")

        st.markdown("**Read Noise (DN)**")
        st.latex(r"\sigma_{\mathrm{read,DN}} = \sqrt{b}")

    with col2:
        st.markdown("**Read Noise (electrons)**")
        st.latex(r"\sigma_{\mathrm{read},e^-} = \frac{\sigma_{\mathrm{read,DN}}}{K} = \frac{\sqrt{b}}{K}")

        st.markdown("**Full Well Capacity**")
        st.latex(r"\mathrm{FWC} = \frac{\mathrm{ADU_{max}}}{K} \quad [e^-]")

        st.markdown("**Dynamic Range**")
        st.latex(r"\mathrm{DR} = 20 \cdot \log_{10}\!\left(\frac{\mathrm{FWC}}{\sigma_{\mathrm{read},e^-}}\right) \quad [\mathrm{dB}]")

    st.markdown("**Maximum SNR** (shot-noise limited):")
    st.latex(r"\mathrm{SNR}_{\max} = \sqrt{\mathrm{FWC}}")

    # ---- 6. Why Vary Exposure, Not Gain? ----
    st.subheader("6. Why Vary Exposure Time, Not Gain?")
    st.markdown(
        "The PTC requires that every data point shares the **same** system gain $K$ "
        "and read noise $\\sigma_{\\mathrm{read}}$. Only varying the illumination level "
        "(via exposure time) satisfies this."
    )

    st.markdown("### The Correct Approach: Exposure Sweep")
    st.markdown(
        "Hold the analog gain constant and vary exposure time. "
        "Longer exposures collect more photons, sweeping the signal from dark to saturation."
    )
    st.latex(r"\mathrm{Var}_i = K \cdot \bar{S}_i + \sigma^2_{\mathrm{read}}")
    st.markdown(
        "Here $K$ and $\\sigma^2_{\\mathrm{read}}$ are **identical for every point**. "
        "The linear fit directly yields the system gain and read noise."
    )
    st.success(
        "Our ExposureSweep dataset (16 exposure times, 3--36 ms, fixed 9 dB gain) "
        "gives K = 3.13 DN/e-, R^2 = 0.996 -- a textbook PTC result."
    )

    st.markdown("### Why Not Sweep Gain?")
    st.markdown(
        "If you sweep analog gain at a fixed exposure, the system gain $K_i$ and "
        "read noise $\\sigma_{\\mathrm{read},i}$ are **different at every data point**:"
    )
    st.latex(r"\mathrm{Var}_i = K_i \cdot \bar{S}_i + \sigma^2_{\mathrm{read},i}")
    st.markdown(
        "Fitting a single line through these points violates the PTC model -- "
        "the slope is not constant. The resulting 'K' would be a meaningless average "
        "across different amplifier configurations."
    )
    st.info(
        "A gain sweep is still useful for sensor characterization: it reveals how "
        "signal, noise, dark level, and saturation change across the gain range. "
        "See the **Gain Sweep Characterization** page for this analysis."
    )

    # ---- 7. The 255 DN Clamp ----
    st.subheader("7. The 255 DN Clamp Effect")
    st.markdown(
        "This sensor stores 12-bit ADC output (0--4095) as 8-bit data via a hard clamp: "
        "`stored = min(native, 255)`. This is **not** a rescale -- 1 stored DN = 1 native DN."
    )

    st.markdown("The clamp creates three zones in the PTC data:")
    st.markdown(
        "| Zone | Condition | Effect |\n"
        "|------|-----------|--------|\n"
        "| **Zone 1: Clean** | Mean << 255 | No pixels hitting clamp. PTC is valid. |\n"
        "| **Zone 2: Partial clamp** | Mean approaching 255 | Some pixels clipped. Variance suppressed. PTC rolls over. |\n"
        "| **Zone 3: Full clamp** | Mean = 255 | All pixels clipped. Variance = 0. No information. |"
    )

    st.warning(
        "FWC measured through the clamp is LIMITED to 255/K electrons -- "
        "this is the capacity through the 8-bit output window, NOT the true sensor "
        "full well capacity. True FWC requires unclamped 12-bit data."
    )

    st.markdown(
        "The saturation threshold (adjustable in the sidebar) controls which pairs "
        "are excluded from the fit. Points in Zone 2 or 3 are rejected to keep the "
        "linear fit clean."
    )

    # ---- 8. Gain Sweep Analysis Method ----
    st.subheader("8. Gain Sweep Analysis Method")
    st.markdown(
        "While a gain sweep cannot produce a valid multi-point PTC (because K changes at "
        "each point), we can still extract K at **each gain level individually** using the "
        "single-point PTC equation."
    )

    st.markdown("**The read-noise bias problem**")
    st.markdown(
        "The naive estimator K_apparent = variance / mean is biased high because "
        "the variance includes read noise:"
    )
    st.latex(
        r"\mathrm{Var} = K \cdot \bar{S} + K^2 \cdot \sigma_{\mathrm{read},e^-}^2 "
        r"\;\;\Longrightarrow\;\; "
        r"K_{\mathrm{apparent}} = \frac{\mathrm{Var}}{\bar{S}} = K + "
        r"\frac{K^2 \cdot \sigma_{\mathrm{read},e^-}^2}{\bar{S}}"
    )
    st.markdown(
        "At 9 dB this overestimate is ~28%. The bias is worse at higher gains "
        "(larger K) and lower signal levels."
    )

    st.markdown("**Quadratic correction for K_true**")
    st.markdown("Treating the single-point PTC equation as a quadratic in K:")
    st.latex(
        r"K^2 \cdot \sigma_{\mathrm{read},e^-}^2 + K \cdot \bar{S} "
        r"- \mathrm{Var} = 0"
    )
    st.markdown("The positive root gives:")
    st.latex(
        r"K_{\mathrm{true}} = \frac{-\bar{S} + \sqrt{\bar{S}^2 + "
        r"4\,\sigma_{\mathrm{read},e^-}^2 \cdot \mathrm{Var}}}"
        r"{2\,\sigma_{\mathrm{read},e^-}^2}"
    )
    st.markdown(
        "This requires an **anchor** value for the input-referred read noise "
        "(sigma_read_e = 2.78 e-), obtained from the ExposureSweep PTC y-intercept."
    )

    st.markdown("**Exponential gain model**")
    st.markdown("K_true grows exponentially with gain:")
    st.latex(r"K(g) = K_0 \cdot 10^{g\,/\,c}")
    st.markdown(
        "where $g$ is gain in dB and $c$ is the dB-per-decade constant. "
        "Fitting in log space:"
    )
    st.latex(r"\log_{10}(K_{\mathrm{true}}) = a + b \cdot g")
    st.markdown("Derived parameters:")
    st.latex(r"K_0 = 10^a \qquad c = \frac{1}{b} \quad (\mathrm{dB\;per\;decade})")
    st.latex(
        r"\mathrm{dB\;per\;doubling} = \frac{\log_{10}(2)}{b} \qquad "
        r"\mathrm{scale\;factor} = \frac{\mathrm{dB_{doubling}}}{6.02}"
    )

    st.markdown("**Key results:**")
    st.markdown(
        f"- K_0 = 1.19 DN/e- (extrapolated to 0 dB; between EMVA LCG = "
        f"{EMVA_SPEC['LCG']['K_dn_per_e']:.3f} and HCG = "
        f"{EMVA_SPEC['HCG']['K_dn_per_e']:.3f})\n"
        "- dB per decade = 21.25, dB per doubling = 6.4\n"
        "- The gain dB scale is close to standard voltage dB "
        "(6.4 dB to double K vs 6.02 standard).\n"
        "- GAINMODE = 1 at all 24 gain levels, confirming a single CG configuration "
        "across the sweep."
    )

    # ---- References ----
    st.subheader("References")
    st.markdown(
        "1. Janesick, J. R., *Photon Transfer: DN → λ*, SPIE Press (2007)\n"
        "2. EMVA Standard 1288, Release 4.0 — Standard for Characterization of Image Sensors and Cameras\n"
        "3. Bohndiek, S. E. et al., *Characterization and Testing of LAS: A Prototype Large Area Sensor*, "
        "IEEE Trans. Nuclear Science **56**(5), 2938–2946 (2009)\n"
        "4. Harvest Imaging, *PTC Tutorial Series*: [harvestimaging.com](https://harvestimaging.com/blog/?p=1034)\n"
        "5. McLean, I. S., *Electronic Imaging in Astronomy: Detectors and Instrumentation*, 2nd ed., Springer (2008)"
    )


# ---------------------------------------------------------------------------
# Page: Frame Explorer
# ---------------------------------------------------------------------------
def page_frame_explorer(settings: dict):
    st.header("Frame Explorer")

    try:
        fit_dict, pair_dicts = run_analysis(settings)
    except ValueError as e:
        st.error(f"Analysis failed: {e}")
        return

    mask = np.array(fit_dict["fit_mask"])
    all_means = np.array([p["mean_signal"] for p in pair_dicts])
    all_variances_raw = np.array([p["variance"] for p in pair_dicts])
    if settings["apply_quant"]:
        all_variances = apply_quantization_correction(all_variances_raw)
    else:
        all_variances = all_variances_raw

    # Get labels
    if LIVE_MODE:
        sweep_folder = settings["sweep_folder"]
        flat_pairs = collect_flat_pairs(sweep_folder)
        labels = [fp.label for fp in flat_pairs]
    else:
        sweep_data = get_demo_sweep(settings)
        labels = sweep_data["labels"]
        fe_data = sweep_data.get("frame_explorer", [])

    # --- Setting selector ---
    sel_idx = st.selectbox(
        "Select exposure/gain level",
        range(len(labels)),
        format_func=lambda i: labels[i],
    )
    pd_pair = pair_dicts[sel_idx]
    pair_used = mask[sel_idx]

    roi = settings["roi"]
    rx, ry, rw, rh = roi

    if LIVE_MODE:
        fp = flat_pairs[sel_idx]
        # --- Load frames ---
        img_a = load_fits_image(fp.bright_a)
        img_b = load_fits_image(fp.bright_b)
        has_dark = fp.dark_a is not None and fp.dark_b is not None
        if has_dark:
            img_dark_a = load_fits_image(fp.dark_a)

        # --- Frame Previews ---
        st.subheader("Frame Previews")
        n_preview_cols = 3 if has_dark else 2
        preview_cols = st.columns(n_preview_cols)

        def _frame_heatmap(img: np.ndarray, title: str) -> go.Figure:
            fig = go.Figure()
            fig.add_trace(go.Heatmap(
                z=img, colorscale="gray", showscale=True,
                colorbar=dict(title="DN"),
            ))
            fig.add_shape(
                type="rect",
                x0=rx, y0=ry, x1=rx + rw, y1=ry + rh,
                line=dict(color="red", width=2),
            )
            fig.update_layout(
                title=title, xaxis_title="X", yaxis_title="Y",
                yaxis=dict(autorange="reversed"),
                height=350, margin=dict(t=40, b=30, l=30, r=30),
            )
            return fig

        with preview_cols[0]:
            st.plotly_chart(
                _frame_heatmap(img_a, f"Bright A: {fp.bright_a.name}"),
                use_container_width=True,
            )
        with preview_cols[1]:
            st.plotly_chart(
                _frame_heatmap(img_b, f"Bright B: {fp.bright_b.name}"),
                use_container_width=True,
            )
        if has_dark:
            with preview_cols[2]:
                st.plotly_chart(
                    _frame_heatmap(img_dark_a, f"Dark: {fp.dark_a.name}"),
                    use_container_width=True,
                )

        # --- Pixel Statistics ---
        st.subheader("Pixel Statistics")

        def _pixel_stats(img: np.ndarray, label: str) -> dict:
            roi_crop = img[ry: ry + rh, rx: rx + rw]
            return {
                "Frame": label,
                "Full Mean": round(float(img.mean()), 2),
                "Full Std": round(float(img.std()), 2),
                "Full Min": int(img.min()),
                "Full Max": int(img.max()),
                "ROI Mean": round(float(roi_crop.mean()), 2),
                "ROI Std": round(float(roi_crop.std()), 2),
                "ROI Min": int(roi_crop.min()),
                "ROI Max": int(roi_crop.max()),
            }

        stats_rows = [
            _pixel_stats(img_a, "Bright A"),
            _pixel_stats(img_b, "Bright B"),
        ]
        if has_dark:
            stats_rows.append(_pixel_stats(img_dark_a, "Dark"))
        st.dataframe(pd.DataFrame(stats_rows), use_container_width=True, hide_index=True)

        # --- Intensity Histogram ---
        st.subheader("ROI Intensity Histogram")
        hist_frame = st.radio(
            "Histogram frame", ["Bright A", "Bright B"] + (["Dark"] if has_dark else []),
            horizontal=True,
        )
        if hist_frame == "Bright A":
            hist_data = img_a[ry: ry + rh, rx: rx + rw].ravel()
        elif hist_frame == "Bright B":
            hist_data = img_b[ry: ry + rh, rx: rx + rw].ravel()
        else:
            hist_data = img_dark_a[ry: ry + rh, rx: rx + rw].ravel()

        fig_hist = go.Figure()
        fig_hist.add_trace(go.Histogram(
            x=hist_data, nbinsx=256, marker_color="steelblue",
            name="Pixel counts",
        ))
        fig_hist.add_vline(x=0, line_dash="dot", line_color="blue",
                           annotation_text="Floor (0)")
        fig_hist.add_vline(x=ADU_MAX, line_dash="dot", line_color="red",
                           annotation_text="Clamp (255)")
        fig_hist.update_layout(
            title=f"ROI Intensity Distribution ({hist_frame})",
            xaxis_title="Pixel Value (DN)", yaxis_title="Count",
            height=350,
        )
        st.plotly_chart(fig_hist, use_container_width=True)

        # --- Difference Image (A - B) ---
        st.subheader("Difference Image (A - B)")
        diff_img = img_a.astype(np.float64) - img_b.astype(np.float64)
        diff_roi = diff_img[ry: ry + rh, rx: rx + rw]

        fig_diff = go.Figure()
        fig_diff.add_trace(go.Heatmap(
            z=diff_roi, colorscale="RdBu_r", zmid=0,
            showscale=True, colorbar=dict(title="DN"),
        ))
        fig_diff.update_layout(
            title=f"A - B (ROI) | mean={diff_roi.mean():.3f}, std={diff_roi.std():.3f}",
            xaxis_title="X (ROI)", yaxis_title="Y (ROI)",
            yaxis=dict(autorange="reversed"),
            height=400,
        )
        st.plotly_chart(fig_diff, use_container_width=True)

    else:
        # DEMO MODE: use precomputed data
        if sel_idx < len(fe_data):
            fe = fe_data[sel_idx]

            st.subheader("Frame Previews")
            st.info(
                "Full-resolution frame previews require local FITS data. "
                "Showing precomputed statistics below."
            )

            # --- Pixel Statistics from demo data ---
            st.subheader("Pixel Statistics")
            stats_rows = [
                {"Frame": f"Bright A ({fe['bright_a_name']})", **fe["stats_a"]},
                {"Frame": f"Bright B ({fe['bright_b_name']})", **fe["stats_b"]},
            ]
            if fe.get("dark_stats"):
                stats_rows.append(
                    {"Frame": f"Dark ({fe['dark_name']})", **fe["dark_stats"]}
                )
            st.dataframe(pd.DataFrame(stats_rows), use_container_width=True, hide_index=True)

            # --- Histogram from precomputed data ---
            st.subheader("ROI Intensity Histogram (Bright A)")
            hist_counts = fe["hist_counts"]
            hist_edges = fe["hist_edges"]
            # Use bar chart with precomputed histogram
            bin_centers = [(hist_edges[i] + hist_edges[i + 1]) / 2 for i in range(len(hist_counts))]
            fig_hist = go.Figure()
            fig_hist.add_trace(go.Bar(
                x=bin_centers, y=hist_counts,
                marker_color="steelblue", name="Pixel counts",
                width=1.0,
            ))
            fig_hist.add_vline(x=0, line_dash="dot", line_color="blue",
                               annotation_text="Floor (0)")
            fig_hist.add_vline(x=ADU_MAX, line_dash="dot", line_color="red",
                               annotation_text="Clamp (255)")
            fig_hist.update_layout(
                title=f"ROI Intensity Distribution (Bright A, {fe['label']})",
                xaxis_title="Pixel Value (DN)", yaxis_title="Count",
                height=350,
            )
            st.plotly_chart(fig_hist, use_container_width=True)

            # --- Diff stats (no image in demo) ---
            st.subheader("Difference Image (A - B) Statistics")
            st.markdown(
                f"**ROI difference:** mean = {fe['diff_roi_mean']:.3f} DN, "
                f"std = {fe['diff_roi_std']:.3f} DN"
            )
        else:
            st.warning("No precomputed frame data available for this setting.")

    # --- PTC Pair Stats Card (works in both modes) ---
    st.subheader("PTC Pair Statistics")
    var_display = all_variances[sel_idx]
    col_s1, col_s2, col_s3, col_s4 = st.columns(4)
    col_s1.metric("Mean Signal", f"{pd_pair['mean_signal']:.2f} DN")
    col_s2.metric("Variance", f"{var_display:.2f} DN^2")
    col_s3.metric("Sat Hi", f"{pd_pair['sat_hi'] * 100:.2f}%")
    col_s4.metric("Sat Lo", f"{pd_pair['sat_lo'] * 100:.2f}%")

    if pair_used:
        st.success("This pair was USED in the PTC linear fit.")
    else:
        st.error("This pair was EXCLUDED from the PTC fit (saturation threshold exceeded).")

    # --- Mini PTC Plot with highlighted point ---
    st.subheader("Position on PTC Curve")
    K = fit_dict["K_slope"]
    intercept = fit_dict["intercept"]

    fig_mini = go.Figure()
    fig_mini.add_trace(go.Scatter(
        x=all_means[mask].tolist(), y=all_variances[mask].tolist(),
        mode="markers", name="Used",
        marker=dict(color="lightblue", size=6),
        hoverinfo="skip",
    ))
    if (~mask).any():
        fig_mini.add_trace(go.Scatter(
            x=all_means[~mask].tolist(), y=all_variances[~mask].tolist(),
            mode="markers", name="Excluded",
            marker=dict(color="lightcoral", size=6, symbol="x"),
            hoverinfo="skip",
        ))
    # Fit line
    used_max = all_means[mask].max() if mask.any() else ADU_MAX
    fit_x = np.linspace(0, min(used_max * 1.2, ADU_MAX), 100)
    fit_y = K * fit_x + intercept
    fig_mini.add_trace(go.Scatter(
        x=fit_x.tolist(), y=fit_y.tolist(),
        mode="lines", name="Fit", line=dict(color="orange", dash="dash"),
    ))
    # Highlight current pair
    fig_mini.add_trace(go.Scatter(
        x=[float(all_means[sel_idx])],
        y=[float(all_variances[sel_idx])],
        mode="markers",
        name=f"Current: {labels[sel_idx]}",
        marker=dict(color="gold", size=14, symbol="star",
                    line=dict(color="black", width=1.5)),
    ))
    fig_mini.add_vline(x=ADU_MAX, line_dash="dot", line_color="gray")
    fig_mini.update_layout(
        title=f"PTC (highlighted: {labels[sel_idx]})",
        xaxis_title="Mean Signal (DN)", yaxis_title="Variance (DN^2)",
        height=380, showlegend=True,
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
    )
    st.plotly_chart(fig_mini, use_container_width=True)

    # --- Bias Comparison (local vs global) -- live mode only ---
    if LIVE_MODE and not settings["is_gain_sweep"]:
        fp = flat_pairs[sel_idx]
        has_dark = fp.dark_a is not None and fp.dark_b is not None
        if has_dark and GLOBAL_BIAS_PATHS[0].exists():
            st.subheader("Bias Comparison: Local Dark vs Global Bias")
            local_bias = make_master_bias([fp.dark_a, fp.dark_b])
            pr_local = mean_var_from_pair(fp.bright_a, fp.bright_b, local_bias, roi)

            global_bias = make_master_bias(GLOBAL_BIAS_PATHS)
            pr_global = mean_var_from_pair(fp.bright_a, fp.bright_b, global_bias, roi)

            bias_rows = [
                {
                    "Bias Method": "Local dark pair",
                    "Mean Signal (DN)": round(pr_local.mean_signal, 2),
                    "Variance (DN^2)": round(pr_local.variance, 2),
                    "Sat Hi (%)": round(pr_local.sat_hi * 100, 2),
                    "Sat Lo (%)": round(pr_local.sat_lo * 100, 2),
                },
                {
                    "Bias Method": "Global master bias",
                    "Mean Signal (DN)": round(pr_global.mean_signal, 2),
                    "Variance (DN^2)": round(pr_global.variance, 2),
                    "Sat Hi (%)": round(pr_global.sat_hi * 100, 2),
                    "Sat Lo (%)": round(pr_global.sat_lo * 100, 2),
                },
            ]
            st.dataframe(pd.DataFrame(bias_rows), use_container_width=True, hide_index=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    st.title("Photon Transfer Curve Analysis")
    st.caption(
        "AstroTracker CMOS sensor characterization | "
        "All analysis in native DN (8-bit clamped at 255, no rescaling)"
    )

    if not LIVE_MODE:
        if not DEMO_DATA_PATH.exists():
            st.error(
                "FITS data files and demo_data.json not found. "
                "Please run with local FITS data or ensure demo_data.json is present."
            )
            st.stop()
        st.sidebar.success("Demo mode: precomputed results")

    settings = build_sidebar()

    # Page navigation
    page = st.sidebar.radio(
        "Page",
        ["PTC Theory & Derivation", "Sensor Overview", "PTC Analysis",
         "Sensor Metrics", "Gain Sweep Characterization", "Frame Explorer"],
        label_visibility="collapsed",
    )

    if page == "PTC Theory & Derivation":
        page_ptc_derivation(settings)
    elif page == "Sensor Overview":
        page_sensor_overview(settings)
    elif page == "PTC Analysis":
        page_ptc_analysis(settings)
    elif page == "Sensor Metrics":
        page_metrics_dashboard(settings)
    elif page == "Gain Sweep Characterization":
        page_gain_sweep(settings)
    elif page == "Frame Explorer":
        page_frame_explorer(settings)


if __name__ == "__main__":
    main()
