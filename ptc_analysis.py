"""
Photon Transfer Curve (PTC) analysis engine.

Processes 8-bit clamped FITS data from an AstroTracker CMOS star tracker.
All analysis is in the stored DN domain (1 stored DN = 1 native 12-bit DN,
clamped at 255). No rescaling is applied.

References:
    - Janesick, "Photon Transfer: DN → λ" (SPIE Press, 2007)
    - EMVA Standard 1288 (v4.0)
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy import stats


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
ADU_MAX = 255                 # 8-bit clamp ceiling
QUANT_VAR = 1.0 / 12.0       # Quantization variance for 1-DN ADC step (DN²)
DEFAULT_ROI = (932, 676, 200, 200)   # (x, y, w, h) — center of 2064×1552
DEFAULT_SAT_THRESH = 0.005    # 0.5 % of ROI pixels at clamp ceiling or floor


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------
@dataclass
class FlatPairInfo:
    """Metadata for a flat-field pair and its associated local dark pair."""
    bright_a: Path
    bright_b: Path
    dark_a: Path | None
    dark_b: Path | None
    label: str                # subfolder name (e.g. "3ms", "90db")


@dataclass
class PairResult:
    """Statistics from a single flat-field pair."""
    path_a: Path
    path_b: Path
    mean_signal: float        # Mean of (A+B)/2 - bias, over ROI (DN)
    variance: float           # 0.5 * var(A-B) over ROI, bias-subtracted (DN^2)
    sat_hi: float             # Fraction of ROI pixels at ADU_MAX in raw frame
    sat_lo: float             # Fraction of ROI pixels at 0 in raw frame
    label: str = ""           # e.g. exposure time or gain label from folder name


@dataclass
class PTCFitResult:
    """Results of the PTC linear fit and derived metrics."""
    # Fit parameters (Var = K_slope * Mean + intercept)
    K_slope: float            # DN/e⁻  (EMVA convention slope)
    intercept: float          # DN² (y-intercept = read noise variance)
    r_squared: float          # Goodness of fit

    # Derived metrics
    conversion_gain: float    # e⁻/DN  (= 1/K_slope, Janesick convention)
    read_noise_dn: float      # DN rms
    read_noise_e: float       # e⁻ rms
    fwc_e: float              # Full well capacity in e⁻ (clamp-limited)
    dynamic_range_db: float   # 20*log10(FWC_e / read_noise_e)
    snr_max: float            # sqrt(FWC_e)

    # Mask info
    n_points_total: int
    n_points_used: int
    fit_mask: np.ndarray = field(repr=False)  # bool array — True = used in fit


# ---------------------------------------------------------------------------
# 1. load_fits_image
# ---------------------------------------------------------------------------
def load_fits_image(path: str | Path) -> np.ndarray:
    """Read a FITS file and return the primary HDU data as float64.

    Parameters
    ----------
    path : str or Path
        Path to the FITS file (uncompressed).

    Returns
    -------
    np.ndarray
        2-D image array (float64).
    """
    with fits.open(str(path)) as hdu_list:
        data = hdu_list[0].data
    if data is None:
        raise ValueError(f"No image data in primary HDU of {path}")
    return data.astype(np.float64)


# ---------------------------------------------------------------------------
# 2. make_master_bias
# ---------------------------------------------------------------------------
def make_master_bias(bias_paths: list[str | Path]) -> np.ndarray:
    """Create a master bias frame by median-combining dark frames.

    Parameters
    ----------
    bias_paths : list of str or Path
        Paths to bias/dark FITS files.

    Returns
    -------
    np.ndarray
        2-D master bias (float64), pixel-wise median.
    """
    if not bias_paths:
        raise ValueError("No bias paths provided")
    stack = np.stack([load_fits_image(p) for p in bias_paths], axis=0)
    return np.median(stack, axis=0)


# ---------------------------------------------------------------------------
# 3. collect_flat_pairs
# ---------------------------------------------------------------------------
def collect_flat_pairs(
    folder: str | Path,
) -> list[FlatPairInfo]:
    """Pair up flat-field FITS files within each subfolder.

    Each subfolder (e.g. "3ms", "90db") contains uncompressed FITS files:
    bright (illuminated) frames and dark frames prefixed with ``dark_``.
    Returns both the bright pair (for PTC analysis) and the dark pair
    (for use as local bias).

    Parameters
    ----------
    folder : str or Path
        Top-level sweep directory (e.g. ExposureSweep_Gain90db/).

    Returns
    -------
    list of FlatPairInfo
        One entry per subfolder with bright pair, dark pair, and label.
    """
    folder = Path(folder)
    pairs: list[FlatPairInfo] = []

    # Natural sort key: extract numeric part from folder names like "3ms", "90db"
    def _sort_key(p: Path) -> float:
        m = re.match(r"(\d+)", p.name)
        return float(m.group(1)) if m else float("inf")

    subdirs = sorted(
        [d for d in folder.iterdir() if d.is_dir()],
        key=_sort_key,
    )

    for subdir in subdirs:
        # Separate bright and dark by filename prefix
        bright = sorted(
            f for f in subdir.glob("*.fits")
            if "_compressed" not in f.name and not f.name.startswith("dark_")
        )
        dark = sorted(
            f for f in subdir.glob("dark_*.fits")
            if "_compressed" not in f.name
        )

        if len(bright) < 2:
            continue  # skip empty or single-frame dirs

        pairs.append(FlatPairInfo(
            bright_a=bright[0],
            bright_b=bright[1],
            dark_a=dark[0] if len(dark) >= 1 else None,
            dark_b=dark[1] if len(dark) >= 2 else None,
            label=subdir.name,
        ))

    return pairs


# ---------------------------------------------------------------------------
# 4. mean_var_from_pair
# ---------------------------------------------------------------------------
def mean_var_from_pair(
    path_a: str | Path,
    path_b: str | Path,
    master_bias: np.ndarray,
    roi: tuple[int, int, int, int] = DEFAULT_ROI,
    adu_max: int = ADU_MAX,
) -> PairResult:
    """Compute mean signal and pair-wise variance from two flat-field frames.

    Uses the Janesick pair method:
        mean  = mean((A + B) / 2) − mean(bias)   over the ROI
        var   = 0.5 * var(A − B)                  over the ROI

    The factor of 0.5 is built in (already halved).

    Saturation fractions are computed on the *raw* (pre-bias-subtraction)
    ROI pixels.

    Parameters
    ----------
    path_a, path_b : str or Path
        Flat-field pair.
    master_bias : np.ndarray
        Master bias frame (float64).
    roi : (x, y, w, h)
        Region of interest in pixel coordinates.
    adu_max : int
        Clamp ceiling (255 for 8-bit data).

    Returns
    -------
    PairResult
    """
    img_a = load_fits_image(path_a)
    img_b = load_fits_image(path_b)

    x, y, w, h = roi
    roi_a = img_a[y : y + h, x : x + w]
    roi_b = img_b[y : y + h, x : x + w]
    roi_bias = master_bias[y : y + h, x : x + w]

    # Saturation fractions on raw data (before bias subtraction)
    n_pixels = roi_a.size
    sat_hi = ((roi_a >= adu_max).sum() + (roi_b >= adu_max).sum()) / (2.0 * n_pixels)
    sat_lo = ((roi_a <= 0).sum() + (roi_b <= 0).sum()) / (2.0 * n_pixels)

    # Bias-subtract
    a_sub = roi_a - roi_bias
    b_sub = roi_b - roi_bias

    # Mean signal: average of both frames
    mean_signal = 0.5 * (a_sub.mean() + b_sub.mean())

    # Pair-wise variance: 0.5 * var(A − B)
    # This cancels fixed-pattern noise (PRNU, DSNU).
    diff = a_sub - b_sub
    variance = 0.5 * np.var(diff)

    return PairResult(
        path_a=Path(path_a),
        path_b=Path(path_b),
        mean_signal=mean_signal,
        variance=variance,
        sat_hi=sat_hi,
        sat_lo=sat_lo,
    )


# ---------------------------------------------------------------------------
# 5. apply_quantization_correction
# ---------------------------------------------------------------------------
def apply_quantization_correction(variance: float | np.ndarray) -> float | np.ndarray:
    """Subtract quantization noise floor from variance.

    The 12-bit ADC has a 1-DN step size (= 1 stored DN because of the clamp),
    giving a quantization variance of 1/12 DN².

    Parameters
    ----------
    variance : float or array
        Raw temporal variance in DN².

    Returns
    -------
    float or array
        Corrected variance, floored at 0.
    """
    return np.maximum(variance - QUANT_VAR, 0.0)


# ---------------------------------------------------------------------------
# 6. fit_ptc
# ---------------------------------------------------------------------------
def fit_ptc(
    means: np.ndarray,
    variances: np.ndarray,
    sat_hi: np.ndarray | None = None,
    sat_lo: np.ndarray | None = None,
    sat_thresh: float = DEFAULT_SAT_THRESH,
) -> PTCFitResult:
    """Fit the linear region of the PTC: Var = K * Mean + σ²_read.

    Saturated pairs (where >sat_thresh fraction of ROI pixels hit the clamp
    ceiling or floor) are excluded from the fit.

    Parameters
    ----------
    means : 1-D array
        Mean signal values (DN) for each pair.
    variances : 1-D array
        Temporal variance values (DN²) for each pair.
    sat_hi, sat_lo : 1-D arrays, optional
        Saturation fractions per pair. If None, all points are used.
    sat_thresh : float
        Maximum allowed saturation fraction for a point to be included.

    Returns
    -------
    PTCFitResult
    """
    means = np.asarray(means, dtype=np.float64)
    variances = np.asarray(variances, dtype=np.float64)

    # Build saturation mask
    n = len(means)
    mask = np.ones(n, dtype=bool)
    if sat_hi is not None:
        mask &= np.asarray(sat_hi) < sat_thresh
    if sat_lo is not None:
        mask &= np.asarray(sat_lo) < sat_thresh

    # Also exclude any points with non-positive mean or variance
    mask &= means > 0
    mask &= variances > 0

    n_used = mask.sum()
    if n_used < 2:
        raise ValueError(
            f"Only {n_used} unsaturated points — need at least 2 for a linear fit."
        )

    # Ordinary least-squares: Var = K * Mean + intercept
    slope, intercept, r_value, _, _ = stats.linregress(
        means[mask], variances[mask]
    )
    r_squared = r_value ** 2

    # Clamp intercept to non-negative (physically, read noise variance ≥ 0)
    intercept = max(intercept, 0.0)

    return derive_metrics(
        K_slope=slope,
        intercept=intercept,
        r_squared=r_squared,
        n_points_total=n,
        n_points_used=int(n_used),
        fit_mask=mask,
    )


# ---------------------------------------------------------------------------
# 7. derive_metrics
# ---------------------------------------------------------------------------
def derive_metrics(
    K_slope: float,
    intercept: float,
    r_squared: float = 0.0,
    n_points_total: int = 0,
    n_points_used: int = 0,
    fit_mask: np.ndarray | None = None,
    adu_max: int = ADU_MAX,
) -> PTCFitResult:
    """Compute sensor metrics from PTC fit parameters.

    Parameters
    ----------
    K_slope : float
        Slope of Var vs Mean (DN/e⁻, EMVA convention).
    intercept : float
        y-intercept (DN²) — read noise variance.
    r_squared : float
        R² of the fit.
    n_points_total, n_points_used : int
        Bookkeeping for the fit mask.
    fit_mask : np.ndarray or None
        Boolean mask of which points were used.
    adu_max : int
        Clamp ceiling in DN.

    Returns
    -------
    PTCFitResult
    """
    if K_slope <= 0:
        raise ValueError(f"K_slope must be positive, got {K_slope:.6g}")

    conversion_gain = 1.0 / K_slope       # e⁻/DN
    read_noise_dn = np.sqrt(intercept)     # DN rms
    read_noise_e = read_noise_dn * conversion_gain   # e⁻ rms

    # FWC: clamp-limited — electrons to reach adu_max
    fwc_e = adu_max * conversion_gain

    # Dynamic range
    if read_noise_e > 0:
        dynamic_range_db = 20.0 * np.log10(fwc_e / read_noise_e)
    else:
        dynamic_range_db = float("inf")

    snr_max = np.sqrt(fwc_e)

    if fit_mask is None:
        fit_mask = np.array([], dtype=bool)

    return PTCFitResult(
        K_slope=K_slope,
        intercept=intercept,
        r_squared=r_squared,
        conversion_gain=conversion_gain,
        read_noise_dn=read_noise_dn,
        read_noise_e=read_noise_e,
        fwc_e=fwc_e,
        dynamic_range_db=dynamic_range_db,
        snr_max=snr_max,
        n_points_total=n_points_total,
        n_points_used=n_points_used,
        fit_mask=fit_mask,
    )


# ---------------------------------------------------------------------------
# High-level pipeline
# ---------------------------------------------------------------------------
def run_ptc_analysis(
    sweep_folder: str | Path,
    bias_paths: list[str | Path] | None = None,
    roi: tuple[int, int, int, int] = DEFAULT_ROI,
    sat_thresh: float = DEFAULT_SAT_THRESH,
    apply_quant_corr: bool = True,
    use_local_dark: bool = True,
) -> tuple[PTCFitResult, list[PairResult]]:
    """Run the full PTC pipeline on a sweep dataset.

    Parameters
    ----------
    sweep_folder : path
        Top-level folder with exposure or gain subfolders.
    bias_paths : list of paths, optional
        Global bias/dark frame FITS files. Used as fallback when local
        dark frames are unavailable and use_local_dark is True, or used
        directly when use_local_dark is False.
    roi : (x, y, w, h)
        Region of interest.
    sat_thresh : float
        Saturation threshold for masking.
    apply_quant_corr : bool
        Whether to subtract 1/12 DN^2 quantization floor.
    use_local_dark : bool
        If True (default), use the dark pair from each subfolder as the
        bias for that exposure level.  Falls back to global bias_paths
        if no local dark pair is available.

    Returns
    -------
    (PTCFitResult, list[PairResult])
    """
    # 1. Global master bias (fallback)
    global_bias: np.ndarray | None = None
    if bias_paths:
        global_bias = make_master_bias(bias_paths)

    # 2. Collect pairs (with local darks)
    flat_pairs = collect_flat_pairs(sweep_folder)
    if not flat_pairs:
        raise ValueError(f"No flat-field pairs found in {sweep_folder}")

    # 3. Compute mean & variance for each pair
    results: list[PairResult] = []
    for fp in flat_pairs:
        # Build local bias from the dark pair in this subfolder
        if use_local_dark and fp.dark_a is not None and fp.dark_b is not None:
            local_bias = make_master_bias([fp.dark_a, fp.dark_b])
        elif global_bias is not None:
            local_bias = global_bias
        else:
            raise ValueError(
                f"No bias available for {fp.label}: no local dark pair "
                "and no global bias_paths provided."
            )

        pr = mean_var_from_pair(fp.bright_a, fp.bright_b, local_bias, roi)
        pr.label = fp.label
        results.append(pr)

    # 4. Quantization correction
    means = np.array([r.mean_signal for r in results])
    variances = np.array([r.variance for r in results])
    if apply_quant_corr:
        variances = apply_quantization_correction(variances)

    sat_hi = np.array([r.sat_hi for r in results])
    sat_lo = np.array([r.sat_lo for r in results])

    # 5. Fit PTC
    fit_result = fit_ptc(means, variances, sat_hi, sat_lo, sat_thresh)

    return fit_result, results


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
def main() -> None:
    """Run PTC analysis on the ExposureSweep dataset and print results."""
    base = Path(__file__).parent / "PhotonTransferCurveArtifacts"

    sweep = base / "ExposureSweep_Gain90db"
    bias_paths = [
        sweep / "BlackImage.fits",
        sweep / "Black_Image2.fits",
    ]

    print("=" * 60)
    print("Photon Transfer Curve Analysis")
    print("Dataset: ExposureSweep_Gain90db")
    print(f"ROI: x={DEFAULT_ROI[0]}, y={DEFAULT_ROI[1]}, "
          f"w={DEFAULT_ROI[2]}, h={DEFAULT_ROI[3]}")
    print(f"Saturation threshold: {DEFAULT_SAT_THRESH * 100:.1f}%")
    print("Bias: local dark pairs per subfolder (global bias as fallback)")
    print("=" * 60)

    fit, pairs = run_ptc_analysis(
        sweep, bias_paths, use_local_dark=True,
    )

    # Per-pair summary table
    print(f"\n{'Label':>6s}  {'Mean DN':>8s}  {'Var DN2':>9s}  "
          f"{'Sat Hi%':>7s}  {'Sat Lo%':>7s}  {'Used?':>5s}")
    print("-" * 55)

    means = np.array([r.mean_signal for r in pairs])
    variances = apply_quantization_correction(np.array([r.variance for r in pairs]))

    for i, r in enumerate(pairs):
        used = "YES" if fit.fit_mask[i] else " NO"
        print(f"{r.label:>6s}  {means[i]:8.2f}  {variances[i]:9.2f}  "
              f"{100 * r.sat_hi:7.2f}  {100 * r.sat_lo:7.2f}  {used:>5s}")

    # Summary metrics
    print(f"\n{'=' * 60}")
    print("PTC Fit Results (stored DN domain, clamp at 255)")
    print(f"{'=' * 60}")
    print(f"  K (slope, DN/e-)     : {fit.K_slope:.6f}")
    print(f"  Intercept (DN2)      : {fit.intercept:.4f}")
    print(f"  R-squared            : {fit.r_squared:.6f}")
    print(f"  Points used / total  : {fit.n_points_used} / {fit.n_points_total}")
    print()
    print(f"  Conversion gain      : {fit.conversion_gain:.2f} e-/DN")
    print(f"  Read noise           : {fit.read_noise_dn:.2f} DN  "
          f"= {fit.read_noise_e:.2f} e- rms")
    print(f"  FWC (clamp-limited)  : {fit.fwc_e:.0f} e-  "
          f"** capped at 255 DN, not true sensor FWC")
    print(f"  Dynamic range        : {fit.dynamic_range_db:.1f} dB")
    print(f"  SNR_max              : {fit.snr_max:.1f}")


if __name__ == "__main__":
    main()
