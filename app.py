"""
Photon Transfer Curve (PTC) Analysis — Streamlit Dashboard

Interactive sensor characterization tool for AstroTracker CMOS star tracker.
All analysis in native DN (stored 8-bit = native 12-bit, clamped at 255).
NO rescaling applied — 1 stored DN = 1 native DN.
"""

from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from astropy.io import fits

from ptc_analysis import (
    ADU_MAX,
    DEFAULT_ROI,
    DEFAULT_SAT_THRESH,
    PairResult,
    PTCFitResult,
    apply_quantization_correction,
    collect_flat_pairs,
    load_fits_image,
    make_master_bias,
    mean_var_from_pair,
    run_ptc_analysis,
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
SWEEP_DIRS = {
    "Exposure Sweep (90 dB gain)": ARTIFACTS / "ExposureSweep_Gain90db",
    "Gain Sweep (10 ms exposure)": ARTIFACTS / "GainSweep_10msExpoTime",
}
GLOBAL_BIAS_PATHS = [
    ARTIFACTS / "ExposureSweep_Gain90db" / "BlackImage.fits",
    ARTIFACTS / "ExposureSweep_Gain90db" / "Black_Image2.fits",
]

# Image dimensions (from sensor spec)
IMG_HEIGHT = 1552
IMG_WIDTH = 2064


# ---------------------------------------------------------------------------
# Sidebar controls
# ---------------------------------------------------------------------------
def build_sidebar() -> dict:
    """Render sidebar controls and return current settings."""
    st.sidebar.title("PTC Analysis Controls")

    dataset_name = st.sidebar.selectbox("Dataset", list(SWEEP_DIRS.keys()))
    sweep_folder = SWEEP_DIRS[dataset_name]
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
        "sweep_folder": sweep_folder,
        "is_gain_sweep": is_gain_sweep,
        "roi": roi,
        "sat_thresh": sat_thresh,
        "apply_quant": apply_quant,
        "use_local_dark": use_local_dark,
    }


# ---------------------------------------------------------------------------
# Cached analysis runner
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
    # Serialize for caching (dataclasses with numpy arrays aren't hashable)
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
    """Wrapper that passes settings to the cached analysis function."""
    bias_str = [str(p) for p in GLOBAL_BIAS_PATHS] if GLOBAL_BIAS_PATHS[0].exists() else None
    return cached_ptc_analysis(
        sweep_folder=str(settings["sweep_folder"]),
        bias_paths_str=bias_str,
        roi=settings["roi"],
        sat_thresh=settings["sat_thresh"],
        apply_quant=settings["apply_quant"],
        use_local_dark=settings["use_local_dark"],
    )


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

    # Dataset info
    sweep = settings["sweep_folder"]
    flat_pairs = collect_flat_pairs(sweep)
    n_pairs = len(flat_pairs)

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Dataset Info")
        st.markdown(f"**Dataset:** {settings['dataset_name']}")
        st.markdown(f"**Path:** `{sweep}`")
        st.markdown(f"**Exposure/gain levels:** {n_pairs}")
        st.markdown(f"**Frames per level:** 2 bright + 2 dark")
        st.markdown(f"**Total usable pairs:** {n_pairs}")

        # Labels
        labels = [fp.label for fp in flat_pairs]
        if settings["is_gain_sweep"]:
            st.markdown(f"**Gain levels:** {', '.join(labels)}")
        else:
            st.markdown(f"**Exposure times:** {', '.join(labels)}")

    with col2:
        st.subheader("Sensor Specifications")
        st.markdown("**Sensor:** AstroTracker CMOS star tracker")
        st.markdown("**Native ADC:** 12-bit (0-4095 DN)")
        st.markdown("**Stored format:** 8-bit via CLAMP (0-255 DN)")
        st.markdown("**Image size:** 2064 x 1552 pixels")
        st.markdown("**Key:** 1 stored DN = 1 native DN (clamped, not rescaled)")

        # Read FITS header from first bright frame
        if flat_pairs:
            first_file = flat_pairs[0].bright_a
            try:
                with fits.open(str(first_file)) as hdul:
                    header = hdul[0].header
                    if len(header) > 0:
                        st.markdown("**FITS Header (sample):**")
                        # Show a few key header cards
                        header_items = []
                        for key in header.keys():
                            if key and key not in ("", "COMMENT", "HISTORY"):
                                header_items.append(f"  {key} = {header[key]}")
                        if header_items:
                            st.code("\n".join(header_items[:15]), language="text")
            except Exception:
                pass

    # Sample image with ROI overlay
    st.subheader("Sample Image with ROI Overlay")
    if flat_pairs:
        sample_path = flat_pairs[len(flat_pairs) // 2].bright_a
        img = load_fits_image(sample_path)

        roi = settings["roi"]
        rx, ry, rw, rh = roi

        fig = go.Figure()
        fig.add_trace(go.Heatmap(
            z=img,
            colorscale="gray",
            showscale=True,
            colorbar=dict(title="DN"),
            reversescale=False,
        ))
        # ROI rectangle
        fig.add_shape(
            type="rect",
            x0=rx, y0=ry, x1=rx + rw, y1=ry + rh,
            line=dict(color="red", width=2),
        )
        fig.add_annotation(
            x=rx + rw / 2, y=ry - 20,
            text=f"ROI ({rw}x{rh})",
            showarrow=False,
            font=dict(color="red", size=12),
        )
        fig.update_layout(
            title=f"Sample: {sample_path.name}",
            xaxis_title="X (pixels)",
            yaxis_title="Y (pixels)",
            yaxis=dict(autorange="reversed"),
            width=900,
            height=700,
        )
        st.plotly_chart(fig, use_container_width=True)
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
        x=fit_x,
        y=fit_y,
        mode="lines",
        name=f"Fit: Var = {K:.4f} * Mean + {intercept:.2f}",
        line=dict(color="orange", dash="dash"),
    ))

    # Vertical clamp ceiling line
    fig.add_vline(
        x=ADU_MAX,
        line_dash="dot",
        line_color="gray",
        annotation_text="255 DN clamp",
        annotation_position="top left",
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
# Page: Gain Sweep Analysis
# ---------------------------------------------------------------------------
def page_gain_sweep(settings: dict):
    st.header("Gain Sweep Analysis")

    if not settings["is_gain_sweep"]:
        st.info(
            "Select the **Gain Sweep (10 ms exposure)** dataset in the sidebar "
            "to use this page. The Exposure Sweep dataset varies exposure time "
            "at a single gain, so per-gain analysis is not applicable."
        )
        return

    sweep_folder = settings["sweep_folder"]
    flat_pairs = collect_flat_pairs(sweep_folder)
    if not flat_pairs:
        st.error("No flat-field pairs found in gain sweep dataset.")
        return

    bias_str = [str(p) for p in GLOBAL_BIAS_PATHS] if GLOBAL_BIAS_PATHS[0].exists() else None

    # Analyze each gain level individually
    gain_labels: list[str] = []
    conv_gains: list[float] = []
    read_noises_e: list[float] = []
    read_noises_dn: list[float] = []
    fwc_es: list[float] = []
    k_slopes: list[float] = []
    r_squareds: list[float] = []
    n_used_list: list[int] = []
    n_total_list: list[int] = []
    errors: list[str] = []

    progress = st.progress(0, text="Analyzing gain levels...")
    for idx, fp in enumerate(flat_pairs):
        progress.progress((idx + 1) / len(flat_pairs), text=f"Analyzing {fp.label}...")
        gain_labels.append(fp.label)

        # Build local bias for this gain
        try:
            if fp.dark_a is not None and fp.dark_b is not None:
                local_bias = make_master_bias([fp.dark_a, fp.dark_b])
            elif bias_str:
                local_bias = make_master_bias([Path(p) for p in bias_str])
            else:
                errors.append(f"{fp.label}: no bias available")
                conv_gains.append(np.nan)
                read_noises_e.append(np.nan)
                read_noises_dn.append(np.nan)
                fwc_es.append(np.nan)
                k_slopes.append(np.nan)
                r_squareds.append(np.nan)
                n_used_list.append(0)
                n_total_list.append(1)
                continue

            # For a single gain level, we only have one bright pair.
            # We can still compute mean/var from this single pair, but can't fit a line.
            # For gain sweep, run the full analysis on the entire sweep to get overall fit,
            # then also extract per-gain-level single-pair statistics.
            from ptc_analysis import mean_var_from_pair

            pr = mean_var_from_pair(
                fp.bright_a, fp.bright_b, local_bias, settings["roi"]
            )
            # With a single pair, we can't do a PTC fit. Store the raw values.
            # The per-gain K comes from the overall sweep fit.
            conv_gains.append(np.nan)  # placeholder, filled from sweep fit
            read_noises_e.append(np.nan)
            read_noises_dn.append(np.nan)
            fwc_es.append(np.nan)
            k_slopes.append(np.nan)
            r_squareds.append(np.nan)
            n_used_list.append(1)
            n_total_list.append(1)

        except Exception as exc:
            errors.append(f"{fp.label}: {exc}")
            conv_gains.append(np.nan)
            read_noises_e.append(np.nan)
            read_noises_dn.append(np.nan)
            fwc_es.append(np.nan)
            k_slopes.append(np.nan)
            r_squareds.append(np.nan)
            n_used_list.append(0)
            n_total_list.append(1)

    progress.empty()

    # Run overall sweep analysis to get the single PTC fit across all gain levels
    st.subheader("Overall Gain Sweep PTC")
    try:
        fit_dict, pair_dicts = run_analysis(settings)
    except ValueError as e:
        st.error(f"Overall sweep analysis failed: {e}")
        return

    # Display overall results
    c1, c2, c3 = st.columns(3)
    c1.metric("K (DN/e-)", f"{fit_dict['K_slope']:.4f}")
    c2.metric("R-squared", f"{fit_dict['r_squared']:.4f}")
    c3.metric("Points used", f"{fit_dict['n_points_used']} / {fit_dict['n_points_total']}")

    # Per-gain pair statistics
    st.subheader("Per-Gain Level Statistics")
    mask = np.array(fit_dict["fit_mask"])
    means = [p["mean_signal"] for p in pair_dicts]
    variances_raw = np.array([p["variance"] for p in pair_dicts])
    if settings["apply_quant"]:
        variances = apply_quantization_correction(variances_raw).tolist()
    else:
        variances = variances_raw.tolist()

    # Extract numeric gain values for plotting
    gain_nums = []
    for label in [p["label"] for p in pair_dicts]:
        try:
            gain_nums.append(int(label.replace("db", "").replace("dB", "")))
        except ValueError:
            gain_nums.append(0)

    # Mean signal vs gain plot
    fig_mean = go.Figure()
    fig_mean.add_trace(go.Scatter(
        x=gain_nums,
        y=means,
        mode="markers+lines",
        name="Mean Signal",
        marker=dict(size=6),
    ))
    fig_mean.add_hline(y=ADU_MAX, line_dash="dot", line_color="gray",
                       annotation_text="255 DN clamp")
    fig_mean.update_layout(
        title="Mean Signal vs Analog Gain",
        xaxis_title="Analog Gain (dB)",
        yaxis_title="Mean Signal (DN)",
        height=400,
    )

    # Variance vs gain plot
    fig_var = go.Figure()
    colors = ["dodgerblue" if m else "red" for m in mask]
    fig_var.add_trace(go.Scatter(
        x=gain_nums,
        y=variances,
        mode="markers+lines",
        name="Variance",
        marker=dict(size=6, color=colors),
    ))
    fig_var.update_layout(
        title="Variance vs Analog Gain",
        xaxis_title="Analog Gain (dB)",
        yaxis_title="Variance (DN^2)",
        height=400,
    )

    col_a, col_b = st.columns(2)
    col_a.plotly_chart(fig_mean, use_container_width=True)
    col_b.plotly_chart(fig_var, use_container_width=True)

    # Saturation fraction vs gain
    fig_sat = go.Figure()
    sat_his = [p["sat_hi"] * 100 for p in pair_dicts]
    sat_los = [p["sat_lo"] * 100 for p in pair_dicts]
    fig_sat.add_trace(go.Scatter(
        x=gain_nums, y=sat_his,
        mode="markers+lines", name="Sat Hi (%)",
        marker=dict(color="red"),
    ))
    fig_sat.add_trace(go.Scatter(
        x=gain_nums, y=sat_los,
        mode="markers+lines", name="Sat Lo (%)",
        marker=dict(color="blue"),
    ))
    fig_sat.add_hline(y=settings["sat_thresh"] * 100, line_dash="dash",
                      line_color="orange", annotation_text="Threshold")
    fig_sat.update_layout(
        title="Saturation Fraction vs Analog Gain",
        xaxis_title="Analog Gain (dB)",
        yaxis_title="Saturation (%)",
        height=400,
    )
    st.plotly_chart(fig_sat, use_container_width=True)

    # PTC plot colored by gain
    st.subheader("PTC Plot (all gain levels)")
    fig_ptc = go.Figure()
    fig_ptc.add_trace(go.Scatter(
        x=[means[i] for i in range(len(means)) if not mask[i]],
        y=[variances[i] for i in range(len(variances)) if not mask[i]],
        mode="markers",
        name="Excluded",
        marker=dict(color="red", size=8, symbol="x"),
        text=[pair_dicts[i]["label"] for i in range(len(pair_dicts)) if not mask[i]],
        hovertemplate="<b>%{text}</b><br>Mean: %{x:.1f} DN<br>Var: %{y:.2f} DN^2<extra></extra>",
    ))
    fig_ptc.add_trace(go.Scatter(
        x=[means[i] for i in range(len(means)) if mask[i]],
        y=[variances[i] for i in range(len(variances)) if mask[i]],
        mode="markers",
        name="Used in fit",
        marker=dict(color="dodgerblue", size=8),
        text=[pair_dicts[i]["label"] for i in range(len(pair_dicts)) if mask[i]],
        hovertemplate="<b>%{text}</b><br>Mean: %{x:.1f} DN<br>Var: %{y:.2f} DN^2<extra></extra>",
    ))
    # Fit line
    K = fit_dict["K_slope"]
    intercept = fit_dict["intercept"]
    used_means = [means[i] for i in range(len(means)) if mask[i]]
    if used_means:
        fit_x = np.linspace(0, min(max(used_means) * 1.2, ADU_MAX), 200)
        fit_y = K * fit_x + intercept
        fig_ptc.add_trace(go.Scatter(
            x=fit_x, y=fit_y,
            mode="lines",
            name=f"Fit: Var = {K:.4f} * Mean + {intercept:.2f}",
            line=dict(color="orange", dash="dash"),
        ))
    fig_ptc.add_vline(x=ADU_MAX, line_dash="dot", line_color="gray",
                      annotation_text="255 DN clamp")
    fig_ptc.update_layout(
        title="Gain Sweep PTC",
        xaxis_title="Mean Signal (DN)",
        yaxis_title="Variance (DN^2)",
        height=550,
    )
    st.plotly_chart(fig_ptc, use_container_width=True)

    # Per-pair table
    st.subheader("Per-Gain Summary Table")
    df = pairs_to_dataframe(pair_dicts, fit_dict, settings["apply_quant"])
    st.dataframe(df, use_container_width=True, hide_index=True)

    if errors:
        st.subheader("Errors")
        for err in errors:
            st.error(err)


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
        "Light incident on a pixel generates photoelectrons. The analog-to-digital "
        "converter (ADC) converts these electrons into Digital Numbers (DN). The "
        "system gain **K** relates the two domains:"
    )
    st.latex(r"S_{\mathrm{DN}} = K \cdot S_{e^-}")
    st.markdown("where:")
    st.markdown(
        "- $S_{\\mathrm{DN}}$ = signal in DN (digital numbers)\n"
        "- $S_{e^-}$ = signal in electrons\n"
        "- $K$ = system gain (DN/e$^-$), the PTC slope\n"
        "- Conversion gain $g = 1/K$ (e$^-$/DN), the Janesick convention"
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

    # ---- 6. Exposure Sweep vs Gain Sweep ----
    st.subheader("6. Exposure Sweep vs Gain Sweep")
    st.markdown(
        "There are two fundamentally different ways to acquire PTC data. "
        "They have **different physical meanings** and **different mathematical validity**."
    )

    # Exposure Sweep
    st.markdown("---")
    st.markdown("### Exposure Sweep (Fixed Gain, Varying Exposure Time)")
    st.markdown(
        "Hold the analog gain constant and vary the exposure time. "
        "Longer exposures collect more photons, increasing the signal level."
    )
    st.markdown(
        "**What varies:** Photon count per pixel (more photons at longer exposure)\n\n"
        "**What is constant:** System gain $K$, read noise $\\sigma_{\\mathrm{read}}$, "
        "all electronics settings"
    )
    st.latex(r"\mathrm{Var}_i = K \cdot \bar{S}_i + \sigma^2_{\mathrm{read}}")
    st.markdown(
        "Here $K$ and $\\sigma^2_{\\mathrm{read}}$ are **the same for every data point**. "
        "The index $i$ runs over exposure times. This is a proper linear model with "
        "one slope and one intercept -- the standard PTC."
    )
    st.success(
        "The exposure sweep is the RIGOROUS PTC measurement. "
        "Each data point samples a different signal level at the same gain, "
        "so the linear fit directly and correctly yields K and read noise for that gain setting."
    )
    st.markdown(
        "**Our ExposureSweep dataset:** 16 exposure times (3--36 ms) at fixed 90 dB gain. "
        "Result: K = 3.13 DN/e$^-$, $R^2$ = 0.996, 9 of 16 points used."
    )

    # Gain Sweep
    st.markdown("---")
    st.markdown("### Gain Sweep (Fixed Exposure, Varying Analog Gain)")
    st.markdown(
        "Hold the exposure time constant and vary the analog gain. "
        "The same number of photons hits the sensor, but the electronic amplification changes."
    )
    st.markdown(
        "**What varies:** System gain $K_i$ (different at each gain setting), "
        "read noise $\\sigma_{\\mathrm{read},i}$ (also gain-dependent)\n\n"
        "**What is constant:** Incident photon count, exposure time"
    )
    st.latex(r"\mathrm{Var}_i = K_i \cdot \bar{S}_i + \sigma^2_{\mathrm{read},i}")
    st.warning(
        "In the gain sweep, BOTH K and sigma_read change at every data point. "
        "The single-line PTC model is violated -- each point has a different slope."
    )

    st.markdown("**Why it still looks like a PTC:**")
    st.markdown(
        "Higher gain amplifies both signal and noise, so higher-gain points have "
        "larger mean AND larger variance. The scatter of (mean, variance) pairs traces "
        "out a curve that *resembles* a PTC, but each point lies on a **different** "
        "underlying line with its own $K_i$."
    )

    st.markdown("**Mathematical comparison:**")
    comp_data = {
        "Property": [
            "PTC equation",
            "K (gain)",
            "Read noise",
            "Linear fit validity",
            "Per-point PTC fit",
            "Best for",
        ],
        "Exposure Sweep": [
            "Var = K * Mean + b",
            "Constant (one K for all points)",
            "Constant (one sigma for all points)",
            "Rigorous -- true linear model",
            "N/A (fit uses all points together)",
            "Characterizing sensor at one gain setting",
        ],
        "Gain Sweep": [
            "Var_i = K_i * Mean_i + b_i",
            "Different at each gain level",
            "Different at each gain level",
            "Approximation -- model is misspecified",
            "Impossible (only 1 pair per gain)",
            "Mapping K, read noise, FWC vs gain",
        ],
    }
    st.dataframe(pd.DataFrame(comp_data), use_container_width=True, hide_index=True)

    st.markdown(
        "**Our GainSweep dataset:** 24 gain levels (30--220 dB) at fixed 10 ms exposure. "
        "The overall PTC fit across all gain levels gives an approximate average $K$. "
        "The more useful analysis is plotting per-gain-level statistics (mean, variance, "
        "saturation fraction) to understand how the sensor behaves across its gain range."
    )

    with st.expander("Deeper: Why can't we fit a per-gain PTC in the gain sweep?"):
        st.markdown(
            "A PTC linear fit requires **at least 2 data points at the same gain setting** "
            "to estimate slope and intercept. Our gain sweep dataset has exactly **1 bright "
            "pair per gain level** (2 bright frames + 2 dark frames = 1 mean/variance point). "
            "With only 1 point, the system of equations is underdetermined:\n\n"
            "$$\\mathrm{Var} = K \\cdot \\bar{S} + b \\quad \\text{(2 unknowns, 1 equation)}$$\n\n"
            "To enable per-gain PTC fitting, you would need to acquire multiple exposure times "
            "at each gain setting -- essentially running an exposure sweep for every gain level."
        )

    # ---- 7. The 255 DN Clamp ----
    st.subheader("7. The 255 DN Clamp Effect")
    st.markdown(
        "This sensor stores 12-bit ADC output (0--4095) as 8-bit data via a hard clamp: "
        "`stored = min(native, 255)`. This is **not** a rescale -- 1 stored DN = 1 native DN."
    )

    st.markdown(
        "The clamp creates three zones in the PTC data:"
    )
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

    # ---- References ----
    st.subheader("References")
    st.markdown(
        "1. Janesick, J. R., *Photon Transfer: DN -> Lambda*, SPIE Press (2007)\n"
        "2. EMVA Standard 1288, Release 4.0 -- Standard for Characterization of Image Sensors\n"
        "3. Bohndiek, S. E. et al., *A guide to photon transfer curve analysis* (2008)\n"
        "4. Harvest Imaging, *PTC Tutorial Series*: harvestimaging.com"
    )


# ---------------------------------------------------------------------------
# Page: Frame Explorer
# ---------------------------------------------------------------------------
def page_frame_explorer(settings: dict):
    st.header("Frame Explorer")

    sweep_folder = settings["sweep_folder"]
    flat_pairs = collect_flat_pairs(sweep_folder)
    if not flat_pairs:
        st.warning("No flat-field pairs found in this dataset.")
        return

    # Run overall PTC analysis for context (mini PTC plot, fit inclusion)
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

    # --- Setting selector ---
    labels = [fp.label for fp in flat_pairs]
    sel_idx = st.selectbox(
        "Select exposure/gain level",
        range(len(labels)),
        format_func=lambda i: labels[i],
    )
    fp = flat_pairs[sel_idx]
    pd_pair = pair_dicts[sel_idx]
    pair_used = mask[sel_idx]

    roi = settings["roi"]
    rx, ry, rw, rh = roi

    # --- Load frames ---
    img_a = load_fits_image(fp.bright_a)
    img_b = load_fits_image(fp.bright_b)
    has_dark = fp.dark_a is not None and fp.dark_b is not None
    if has_dark:
        img_dark_a = load_fits_image(fp.dark_a)

    # --- Frame Previews (Bright A, Bright B, Dark) ---
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
        roi_crop = img[ry : ry + rh, rx : rx + rw]
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
        hist_data = img_a[ry : ry + rh, rx : rx + rw].ravel()
    elif hist_frame == "Bright B":
        hist_data = img_b[ry : ry + rh, rx : rx + rw].ravel()
    else:
        hist_data = img_dark_a[ry : ry + rh, rx : rx + rw].ravel()

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
    diff_roi = diff_img[ry : ry + rh, rx : rx + rw]

    fig_diff = go.Figure()
    fig_diff.add_trace(go.Heatmap(
        z=diff_roi,
        colorscale="RdBu_r",
        zmid=0,
        showscale=True,
        colorbar=dict(title="DN"),
    ))
    fig_diff.update_layout(
        title=f"A - B (ROI) | mean={diff_roi.mean():.3f}, std={diff_roi.std():.3f}",
        xaxis_title="X (ROI)", yaxis_title="Y (ROI)",
        yaxis=dict(autorange="reversed"),
        height=400,
    )
    st.plotly_chart(fig_diff, use_container_width=True)

    # --- PTC Pair Stats Card ---
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
    # All points (faded)
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
        name=f"Current: {fp.label}",
        marker=dict(color="gold", size=14, symbol="star",
                    line=dict(color="black", width=1.5)),
    ))
    fig_mini.add_vline(x=ADU_MAX, line_dash="dot", line_color="gray")
    fig_mini.update_layout(
        title=f"PTC (highlighted: {fp.label})",
        xaxis_title="Mean Signal (DN)", yaxis_title="Variance (DN^2)",
        height=380, showlegend=True,
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
    )
    st.plotly_chart(fig_mini, use_container_width=True)

    # --- Bias Comparison (local vs global) ---
    if not settings["is_gain_sweep"] and has_dark and GLOBAL_BIAS_PATHS[0].exists():
        st.subheader("Bias Comparison: Local Dark vs Global Bias")
        # Local dark bias result (already computed in the current analysis)
        local_bias = make_master_bias([fp.dark_a, fp.dark_b])
        pr_local = mean_var_from_pair(fp.bright_a, fp.bright_b, local_bias, roi)

        # Global bias result
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

    settings = build_sidebar()

    # Page navigation
    page = st.sidebar.radio(
        "Page",
        ["PTC Theory & Derivation", "Sensor Overview", "PTC Analysis",
         "Sensor Metrics", "Gain Sweep Analysis", "Frame Explorer"],
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
    elif page == "Gain Sweep Analysis":
        page_gain_sweep(settings)
    elif page == "Frame Explorer":
        page_frame_explorer(settings)


if __name__ == "__main__":
    main()
