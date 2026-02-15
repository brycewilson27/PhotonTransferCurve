"""Generate precomputed demo data for Streamlit Cloud deployment."""
import json
import numpy as np
from pathlib import Path
from ptc_analysis import (
    collect_flat_pairs, run_ptc_analysis, apply_quantization_correction,
    load_fits_image, make_master_bias, mean_var_from_pair,
    analyze_gain_sweep_per_point, fit_gain_model,
    DEFAULT_ROI, DEFAULT_SAT_THRESH, ADU_MAX
)
from astropy.io import fits as afits

base = Path("PhotonTransferCurveArtifacts")
roi = DEFAULT_ROI
sat_thresh = DEFAULT_SAT_THRESH
demo = {}


def process_sweep(sweep_path, bias_paths, sweep_key, is_gain_sweep=False):
    """Process one sweep dataset and return demo dict."""
    print(f"Processing {sweep_key}...")
    fit, pairs = run_ptc_analysis(
        sweep_path, bias_paths, roi, sat_thresh,
        apply_quant_corr=True, use_local_dark=True,
    )

    pair_dicts = []
    for pr in pairs:
        pair_dicts.append({
            "label": pr.label,
            "mean_signal": float(pr.mean_signal),
            "variance": float(pr.variance),
            "sat_hi": float(pr.sat_hi),
            "sat_lo": float(pr.sat_lo),
            "path_a": pr.path_a.name,
            "path_b": pr.path_b.name,
        })

    fit_dict = {
        "K_slope": float(fit.K_slope),
        "intercept": float(fit.intercept),
        "r_squared": float(fit.r_squared),
        "conversion_gain": float(fit.conversion_gain),
        "read_noise_dn": float(fit.read_noise_dn),
        "read_noise_e": float(fit.read_noise_e),
        "fwc_e": float(fit.fwc_e),
        "dynamic_range_db": float(fit.dynamic_range_db),
        "snr_max": float(fit.snr_max),
        "n_points_total": int(fit.n_points_total),
        "n_points_used": int(fit.n_points_used),
        "fit_mask": fit.fit_mask.tolist(),
    }

    flat_pairs = collect_flat_pairs(sweep_path)
    labels = [fp.label for fp in flat_pairs]

    result = {
        "fit": fit_dict,
        "pairs": pair_dicts,
        "labels": labels,
        "n_pairs": len(flat_pairs),
    }

    # Sample image for sensor overview (exposure sweep only)
    if not is_gain_sweep:
        sample_fp = flat_pairs[len(flat_pairs) // 2]
        sample_img = load_fits_image(sample_fp.bright_a)
        ds = 4
        result["sample_image_ds"] = sample_img[::ds, ::ds].tolist()
        result["sample_image_ds_factor"] = ds
        result["sample_image_name"] = sample_fp.bright_a.name

        # FITS header
        with afits.open(str(flat_pairs[0].bright_a)) as hdul:
            hdr = hdul[0].header
            header_items = []
            for key in hdr.keys():
                if key and key not in ("", "COMMENT", "HISTORY"):
                    header_items.append(f"{key} = {hdr[key]}")
            result["header_text"] = "\n".join(header_items[:15])

    # Frame explorer data for each pair
    rx, ry, rw, rh = roi
    frame_explorer = []
    for fp in flat_pairs:
        img_a = load_fits_image(fp.bright_a)
        img_b = load_fits_image(fp.bright_b)
        roi_a = img_a[ry:ry + rh, rx:rx + rw]
        roi_b = img_b[ry:ry + rh, rx:rx + rw]

        has_dark = fp.dark_a is not None and fp.dark_b is not None
        dark_stats = None
        if has_dark:
            img_dark = load_fits_image(fp.dark_a)
            roi_dark = img_dark[ry:ry + rh, rx:rx + rw]
            dark_stats = {
                "Full Mean": round(float(img_dark.mean()), 2),
                "Full Std": round(float(img_dark.std()), 2),
                "Full Min": int(img_dark.min()),
                "Full Max": int(img_dark.max()),
                "ROI Mean": round(float(roi_dark.mean()), 2),
                "ROI Std": round(float(roi_dark.std()), 2),
                "ROI Min": int(roi_dark.min()),
                "ROI Max": int(roi_dark.max()),
            }

        diff = img_a.astype(np.float64) - img_b.astype(np.float64)
        diff_roi = diff[ry:ry + rh, rx:rx + rw]
        hist_counts, hist_edges = np.histogram(roi_a.ravel(), bins=256, range=(0, 256))

        fe = {
            "label": fp.label,
            "bright_a_name": fp.bright_a.name,
            "bright_b_name": fp.bright_b.name,
            "dark_name": fp.dark_a.name if fp.dark_a else None,
            "stats_a": {
                "Full Mean": round(float(img_a.mean()), 2),
                "Full Std": round(float(img_a.std()), 2),
                "Full Min": int(img_a.min()),
                "Full Max": int(img_a.max()),
                "ROI Mean": round(float(roi_a.mean()), 2),
                "ROI Std": round(float(roi_a.std()), 2),
                "ROI Min": int(roi_a.min()),
                "ROI Max": int(roi_a.max()),
            },
            "stats_b": {
                "Full Mean": round(float(img_b.mean()), 2),
                "Full Std": round(float(img_b.std()), 2),
                "Full Min": int(img_b.min()),
                "Full Max": int(img_b.max()),
                "ROI Mean": round(float(roi_b.mean()), 2),
                "ROI Std": round(float(roi_b.std()), 2),
                "ROI Min": int(roi_b.min()),
                "ROI Max": int(roi_b.max()),
            },
            "dark_stats": dark_stats,
            "diff_roi_mean": round(float(diff_roi.mean()), 3),
            "diff_roi_std": round(float(diff_roi.std()), 3),
            "hist_counts": hist_counts.tolist(),
            "hist_edges": hist_edges.tolist(),
            # Downsampled images for preview
            "img_a_ds": img_a[::8, ::8].tolist(),
            "img_b_ds": img_b[::8, ::8].tolist(),
            "diff_roi": diff_roi[::2, ::2].tolist(),
        }
        if has_dark:
            fe["img_dark_ds"] = load_fits_image(fp.dark_a)[::8, ::8].tolist()
        frame_explorer.append(fe)
        print(f"  {fp.label} done")

    result["frame_explorer"] = frame_explorer
    return result


if __name__ == "__main__":
    # Exposure Sweep
    sweep_exp = base / "ExposureSweep_Gain9db"
    bias_paths = [sweep_exp / "BlackImage.fits", sweep_exp / "Black_Image2.fits"]
    demo["exposure_sweep"] = process_sweep(sweep_exp, bias_paths, "ExposureSweep")

    # Gain Sweep
    sweep_gain = base / "GainSweep_10msExpoTime"
    demo["gain_sweep"] = process_sweep(sweep_gain, None, "GainSweep", is_gain_sweep=True)

    # Per-point gain analysis (uses ExposureSweep read_noise_e as anchor)
    sigma_read_e_anchor = demo["exposure_sweep"]["fit"]["read_noise_e"]
    print(f"\nRunning per-point gain analysis (anchor={sigma_read_e_anchor:.2f} e-)...")
    gp_results = analyze_gain_sweep_per_point(
        sweep_gain, roi, sigma_read_e_anchor=sigma_read_e_anchor, apply_quant_corr=True
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
    demo["gain_sweep"]["per_point_analysis"] = per_point
    print(f"  {len(per_point)} gain points processed")

    # Gain model fit
    gm = fit_gain_model(gp_results, sat_thresh=DEFAULT_SAT_THRESH)
    demo["gain_sweep"]["gain_model"] = {
        "K_0": float(gm.K_0),
        "dB_per_decade": float(gm.dB_per_decade),
        "dB_per_doubling": float(gm.dB_per_doubling),
        "scale_factor": float(gm.scale_factor),
        "r_squared": float(gm.r_squared),
        "n_points_used": gm.n_points_used,
        "gain_range": list(gm.gain_range),
        "sigma_read_e_mean": float(gm.sigma_read_e_mean),
    }
    print(f"  Gain model: K_0={gm.K_0:.4f}, R^2={gm.r_squared:.6f}, "
          f"{gm.n_points_used} points")

    # Metadata
    demo["metadata"] = {
        "roi": list(roi),
        "sat_thresh": sat_thresh,
        "adu_max": ADU_MAX,
        "img_height": 1552,
        "img_width": 2064,
        "apply_quant": True,
        "use_local_dark": True,
    }

    with open("demo_data.json", "w") as f:
        json.dump(demo, f)

    import os
    size_mb = os.path.getsize("demo_data.json") / 1048576
    print(f"\ndemo_data.json written: {size_mb:.1f} MB")
    print("Done!")
