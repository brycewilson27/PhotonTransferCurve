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

# Lazy imports: astropy and scipy are only needed for FITS I/O and fitting,
# not for constants/dicts. Defer so demo mode (Streamlit Cloud) works even
# if these packages fail to install.
fits = None   # astropy.io.fits -- loaded on first use
stats = None  # scipy.stats -- loaded on first use


def _ensure_fits():
    """Lazy-load astropy.io.fits."""
    global fits
    if fits is None:
        from astropy.io import fits as _fits
        fits = _fits


def _ensure_stats():
    """Lazy-load scipy.stats."""
    global stats
    if stats is None:
        from scipy import stats as _stats
        stats = _stats


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
ADU_MAX = 255                 # 8-bit clamp ceiling
QUANT_VAR = 1.0 / 12.0       # Quantization variance for 1-DN ADC step (DN²)
DEFAULT_ROI = (932, 676, 200, 200)   # (x, y, w, h) — center of 2064×1552
DEFAULT_SAT_THRESH = 0.005    # 0.5 % of ROI pixels at clamp ceiling or floor


# ---------------------------------------------------------------------------
# EMVA 1288 Reference Specifications (Sony IMX900-AMR-C, FRAMOS report)
# ---------------------------------------------------------------------------
# Two operating points from the FRAMOS EMVA1288-compliant characterization.
# These are the manufacturer/lab-grade specs at 12-bit output, 536 nm LED.
# Our AstroTracker measurements are at 9 dB gain with 8-bit clamped output,
# so direct comparison is not possible -- but these provide ground truth for
# the underlying sensor capability.
#
# Source: FRAMOS Technologies Inc. EMVA1288 datasheets for IMX900-AMR-C
#   - LCG report: Analog Gain 0 dB, Conversion Gain Low
#   - HCG report: Analog Gain 3 dB, Conversion Gain High

EMVA_SPEC = {
    "sensor_model": "Sony IMX900-AMR-C",
    "sensor_type": "Global Shutter CMOS",
    "resolution": "2064(H) x 1552(V)",
    "pixel_size_um": 2.25,
    "diagonal_mm": 5.81,
    "adc_bit_depth": 12,
    "test_illumination": "536 nm LED (32 nm width)",
    "test_temperature_C": 36,
    "test_standard": "EMVA 1288 v3.1",
    "test_lab": "FRAMOS Technologies Inc.",
    "LCG": {
        "label": "Low Conversion Gain (0 dB)",
        "analog_gain_dB": 0,
        "conversion_gain_mode": "Low",
        "K_dn_per_e": 0.399,           # DN/e-  (system gain)
        "conv_gain_e_per_dn": 2.505,    # e-/DN  (1/K)
        "qe_percent": 91.25,            # quantum efficiency at 536 nm
        "read_noise_e": 5.560,           # e- rms (temporal dark noise)
        "read_noise_dn": 2.238,          # DN rms
        "fwc_e": 9458,                   # saturation capacity (electrons)
        "fwc_photons": 10365,            # saturation capacity (photons)
        "snr_max": 97,                   # peak SNR (linear)
        "snr_max_db": 39.76,             # peak SNR (dB)
        "dynamic_range_db": 63.8,        # DR in dB
        "dynamic_range_bits": 10.6,      # DR in bits
        "sensitivity_threshold_photons": 6.692,  # min detectable photons
        "sensitivity_threshold_e": 6.106,         # min detectable electrons
        "dsnu_e": 3.7,                   # dark signal non-uniformity (e-)
        "dsnu_dn": 1.5,                  # DSNU (DN)
        "prnu_percent": 0.7,             # photo-response non-uniformity (%)
        "linearity_error_min_pct": -1.133,
        "linearity_error_max_pct": 0.595,
        "black_level_dn": 240,           # digital black level offset (LSB)
    },
    "HCG": {
        "label": "High Conversion Gain (3 dB)",
        "analog_gain_dB": 3,
        "conversion_gain_mode": "High",
        "K_dn_per_e": 1.673,
        "conv_gain_e_per_dn": 0.598,
        "qe_percent": 88.47,
        "read_noise_e": 1.954,
        "read_noise_dn": 3.281,
        "fwc_e": 2183,
        "fwc_photons": 2468,
        "snr_max": 47,
        "snr_max_db": 33.39,
        "dynamic_range_db": 59.0,
        "dynamic_range_bits": 9.8,
        "sensitivity_threshold_photons": 2.782,
        "sensitivity_threshold_e": 2.461,
        "dsnu_e": 1.3,
        "dsnu_dn": 2.1,
        "prnu_percent": 0.8,
        "linearity_error_min_pct": -0.550,
        "linearity_error_max_pct": 0.294,
        "dark_current_e_per_s": 1.473,
        "dark_current_dn_per_s": 2.464,
        "black_level_dn": 240,
    },
}


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
    label: str                # subfolder name (e.g. "3ms", "9db")


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


@dataclass
class GainPointResult:
    """Per-gain-level analysis from the gain sweep dataset.

    Each instance represents one gain setting with per-point K estimates
    derived from a single bright pair and dark pair.
    """
    gain_db: float            # Analog gain setting (dB)
    mean_signal: float        # Mean signal in ROI after bias subtraction (DN)
    pair_variance: float      # 0.5 * var(A-B) in ROI, bias-subtracted (DN^2)
    dark_pair_variance: float # 0.5 * var(darkA - darkB) in ROI (DN^2)
    K_apparent: float         # variance / mean -- biased high by read noise (DN/e-)
    K_true: float             # Read-noise-corrected K (DN/e-), or NaN if invalid
    read_noise_dn: float      # sqrt(dark_pair_variance) (DN rms)
    read_noise_e: float       # read_noise_dn / K_true (e- rms), or NaN
    n_electrons: float        # mean_signal / K_true (e-), or NaN
    sat_hi: float             # Fraction of ROI pixels at ADU_MAX
    sat_lo: float             # Fraction of ROI pixels at 0
    # FITS header metadata
    gain_mode: int = 0        # GAINMODE from FITS header
    exptime_ns: int = 0       # EXPTIME from FITS header (nanoseconds)


@dataclass
class GainModelResult:
    """Exponential model fit: K_true vs gain_dB.

    Model: log10(K_true) = a + b * gain_dB
    => K_true(g) = K_0 * 10^(g / dB_per_decade)

    The scale_factor quantifies the ratio of the measured dB-per-doubling
    to the standard 6.02 dB (voltage) or 3.01 dB (power).
    """
    K_0: float                # K at 0 dB (extrapolated intercept, DN/e-)
    dB_per_decade: float      # dB for a 10x increase in K
    dB_per_doubling: float    # dB for a 2x increase in K
    scale_factor: float       # dB_per_doubling / 6.02 (ratio to std voltage dB)
    r_squared: float          # R^2 of log-linear fit
    n_points_used: int        # Number of gain levels in the fit
    gain_range: tuple[float, float]  # (min, max) gain_dB used in fit
    sigma_read_e_mean: float  # Mean input-referred read noise across used points (e-)


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
    _ensure_fits()
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

    Each subfolder (e.g. "3ms", "9db") contains uncompressed FITS files:
    bright (illuminated) frames and dark frames prefixed with ``dark_``.
    Returns both the bright pair (for PTC analysis) and the dark pair
    (for use as local bias).

    Parameters
    ----------
    folder : str or Path
        Top-level sweep directory (e.g. ExposureSweep_Gain9db/).

    Returns
    -------
    list of FlatPairInfo
        One entry per subfolder with bright pair, dark pair, and label.
    """
    folder = Path(folder)
    pairs: list[FlatPairInfo] = []

    # Natural sort key: extract numeric part from folder names like "3ms", "9db"
    def _sort_key(p: Path) -> float:
        m = re.match(r"(\d+\.?\d*)", p.name)
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
    _ensure_stats()
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
# Per-point gain sweep analysis
# ---------------------------------------------------------------------------
def analyze_gain_sweep_per_point(
    sweep_folder: str | Path,
    roi: tuple[int, int, int, int] = DEFAULT_ROI,
    sigma_read_e_anchor: float = 2.78,
    apply_quant_corr: bool = True,
) -> list[GainPointResult]:
    """Compute per-gain-level K estimates from the gain sweep dataset.

    For each gain level, computes K_apparent (naive variance/mean) and
    K_true (corrected for read-noise bias via quadratic solution).

    Parameters
    ----------
    sweep_folder : path
        Top-level GainSweep directory (e.g. GainSweep_10msExpoTime/).
    roi : (x, y, w, h)
        Region of interest.
    sigma_read_e_anchor : float
        Input-referred read noise in electrons from the exposure sweep PTC.
        Used for the K_true correction (default 2.78 e- from ExposureSweep).
    apply_quant_corr : bool
        Whether to subtract 1/12 DN^2 quantization floor from variances.

    Returns
    -------
    list[GainPointResult]
        One entry per gain level, sorted by gain_db.
    """
    sweep_folder = Path(sweep_folder)
    flat_pairs = collect_flat_pairs(sweep_folder)
    if not flat_pairs:
        raise ValueError(f"No flat-field pairs found in {sweep_folder}")

    x, y, w, h = roi
    results: list[GainPointResult] = []

    for fp in flat_pairs:
        # --- Local dark bias from subfolder dark pair ---
        if fp.dark_a is None or fp.dark_b is None:
            continue  # skip if no dark pair
        local_bias = make_master_bias([fp.dark_a, fp.dark_b])

        # --- Bright pair: mean & variance ---
        pr = mean_var_from_pair(fp.bright_a, fp.bright_b, local_bias, roi)

        pair_var = pr.variance
        if apply_quant_corr:
            pair_var = apply_quantization_correction(pair_var)

        # --- Dark pair-diff read noise ---
        dark_a_img = load_fits_image(fp.dark_a)
        dark_b_img = load_fits_image(fp.dark_b)
        dark_roi_a = dark_a_img[y: y + h, x: x + w]
        dark_roi_b = dark_b_img[y: y + h, x: x + w]
        dark_pair_var = 0.5 * np.var(dark_roi_a - dark_roi_b)
        read_noise_dn = np.sqrt(dark_pair_var)

        # --- Extract gain_db from folder label ---
        label = fp.label
        m = re.match(r"(\d+\.?\d*)", label)
        gain_db = float(m.group(1)) if m else 0.0

        # --- K_apparent: naive variance / mean (biased by read noise) ---
        mean_sig = pr.mean_signal
        if mean_sig > 0 and pair_var > 0:
            K_apparent = pair_var / mean_sig
        else:
            K_apparent = float("nan")

        # --- K_true: read-noise-corrected via quadratic ---
        # PTC single-point: var = K*mean + K^2 * sigma_read_e^2
        # => K^2 * sigma_read_e^2 + K * mean - var = 0
        # => K = (-mean + sqrt(mean^2 + 4*sigma_read_e^2*var)) / (2*sigma_read_e^2)
        sigma2 = sigma_read_e_anchor ** 2
        discriminant = mean_sig ** 2 + 4.0 * sigma2 * pair_var
        if discriminant > 0 and mean_sig > 0 and pair_var > 0:
            K_true = (-mean_sig + np.sqrt(discriminant)) / (2.0 * sigma2)
            if K_true <= 0:
                K_true = float("nan")
        else:
            K_true = float("nan")

        # --- Derived quantities ---
        if np.isfinite(K_true) and K_true > 0:
            n_electrons = mean_sig / K_true
            read_noise_e = read_noise_dn / K_true
        else:
            n_electrons = float("nan")
            read_noise_e = float("nan")

        # --- FITS header metadata from one bright frame ---
        gain_mode = 0
        exptime_ns = 0
        try:
            _ensure_fits()
            with fits.open(str(fp.bright_a)) as hdul:
                hdr = hdul[0].header
                gain_mode = int(hdr.get("GAINMODE", 0))
                exptime_ns = int(hdr.get("EXPTIME", 0))
        except Exception:
            pass  # header read failure is non-fatal

        results.append(GainPointResult(
            gain_db=gain_db,
            mean_signal=mean_sig,
            pair_variance=pair_var,
            dark_pair_variance=dark_pair_var,
            K_apparent=K_apparent,
            K_true=K_true,
            read_noise_dn=read_noise_dn,
            read_noise_e=read_noise_e,
            n_electrons=n_electrons,
            sat_hi=pr.sat_hi,
            sat_lo=pr.sat_lo,
            gain_mode=gain_mode,
            exptime_ns=exptime_ns,
        ))

    # Sort by gain_db
    results.sort(key=lambda r: r.gain_db)
    return results


# ---------------------------------------------------------------------------
# Gain model fitting
# ---------------------------------------------------------------------------
def fit_gain_model(
    gain_points: list[GainPointResult],
    sat_thresh: float = 0.005,
) -> GainModelResult:
    """Fit an exponential model to K_true vs gain_dB.

    Model: log10(K_true) = a + b * gain_dB
    => K_true(g) = 10^a * 10^(b*g) = K_0 * 10^(g / dB_per_decade)

    Parameters
    ----------
    gain_points : list[GainPointResult]
        Per-gain-level results from analyze_gain_sweep_per_point().
    sat_thresh : float
        Maximum sat_hi to include a point (default 0.5%).

    Returns
    -------
    GainModelResult
    """
    # Filter to unsaturated points with valid K_true
    used = [
        gp for gp in gain_points
        if gp.sat_hi < sat_thresh and np.isfinite(gp.K_true) and gp.K_true > 0
    ]
    if len(used) < 2:
        raise ValueError(
            f"Only {len(used)} valid points after filtering -- need at least 2."
        )

    gain_arr = np.array([gp.gain_db for gp in used], dtype=np.float64)
    log_K = np.log10(np.array([gp.K_true for gp in used], dtype=np.float64))

    # Linear fit: log10(K_true) = a + b * gain_dB
    _ensure_stats()
    slope_b, intercept_a, r_value, _, _ = stats.linregress(gain_arr, log_K)
    r_squared = r_value ** 2

    # Derived parameters
    K_0 = 10.0 ** intercept_a           # K at 0 dB (extrapolated)
    dB_per_decade = 1.0 / slope_b       # dB for 10x K increase
    dB_per_doubling = np.log10(2.0) / slope_b  # dB for 2x K increase
    scale_factor = dB_per_doubling / 6.02  # ratio to standard voltage dB

    # Mean input-referred read noise across used points
    rne_vals = [gp.read_noise_e for gp in used if np.isfinite(gp.read_noise_e)]
    sigma_read_e_mean = float(np.mean(rne_vals)) if rne_vals else float("nan")

    return GainModelResult(
        K_0=K_0,
        dB_per_decade=dB_per_decade,
        dB_per_doubling=dB_per_doubling,
        scale_factor=scale_factor,
        r_squared=r_squared,
        n_points_used=len(used),
        gain_range=(float(gain_arr.min()), float(gain_arr.max())),
        sigma_read_e_mean=sigma_read_e_mean,
    )


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
def main() -> None:
    """Run PTC analysis on ExposureSweep, then per-point gain sweep analysis."""
    base = Path(__file__).parent / "PhotonTransferCurveArtifacts"

    # ---------------------------------------------------------------
    # Part 1: Exposure Sweep PTC (existing)
    # ---------------------------------------------------------------
    sweep = base / "ExposureSweep_Gain9db"
    bias_paths = [
        sweep / "BlackImage.fits",
        sweep / "Black_Image2.fits",
    ]

    print("=" * 70)
    print("Photon Transfer Curve Analysis")
    print("Dataset: ExposureSweep_Gain9db")
    print(f"ROI: x={DEFAULT_ROI[0]}, y={DEFAULT_ROI[1]}, "
          f"w={DEFAULT_ROI[2]}, h={DEFAULT_ROI[3]}")
    print(f"Saturation threshold: {DEFAULT_SAT_THRESH * 100:.1f}%")
    print("Bias: local dark pairs per subfolder (global bias as fallback)")
    print("=" * 70)

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
    print(f"\n{'=' * 70}")
    print("PTC Fit Results (stored DN domain, clamp at 255)")
    print(f"{'=' * 70}")
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

    # ---------------------------------------------------------------
    # Part 2: Per-Point Gain Sweep Analysis (Session A)
    # ---------------------------------------------------------------
    gain_sweep = base / "GainSweep_10msExpoTime"
    sigma_read_e = fit.read_noise_e  # anchor from exposure sweep

    print(f"\n\n{'=' * 70}")
    print("Per-Point Gain Sweep Analysis")
    print(f"Dataset: GainSweep_10msExpoTime")
    print(f"sigma_read_e anchor: {sigma_read_e:.2f} e- (from ExposureSweep PTC)")
    print(f"{'=' * 70}")

    gp_results = analyze_gain_sweep_per_point(
        gain_sweep,
        roi=DEFAULT_ROI,
        sigma_read_e_anchor=sigma_read_e,
        apply_quant_corr=True,
    )

    # Per-gain summary table
    print(f"\n{'Gain':>5s}  {'Mean':>7s}  {'Var':>8s}  {'K_app':>7s}  "
          f"{'K_true':>7s}  {'RN_dn':>6s}  {'RN_e':>6s}  "
          f"{'n_e':>7s}  {'SatHi%':>6s}")
    print("-" * 75)

    for gp in gp_results:
        # Format NaN fields as dashes
        def fmt(v, w, d):
            if np.isfinite(v):
                return f"{v:{w}.{d}f}"
            return "-".rjust(w)

        print(f"{gp.gain_db:>5.1f}  "
              f"{fmt(gp.mean_signal, 7, 2)}  "
              f"{fmt(gp.pair_variance, 8, 2)}  "
              f"{fmt(gp.K_apparent, 7, 3)}  "
              f"{fmt(gp.K_true, 7, 3)}  "
              f"{fmt(gp.read_noise_dn, 6, 2)}  "
              f"{fmt(gp.read_noise_e, 6, 2)}  "
              f"{fmt(gp.n_electrons, 7, 1)}  "
              f"{100 * gp.sat_hi:6.2f}")

    # ---------------------------------------------------------------
    # Verification checks
    # ---------------------------------------------------------------
    print(f"\n{'=' * 70}")
    print("Verification Checks")
    print(f"{'=' * 70}")

    # Check 1: K_true at 9 dB vs exposure sweep K
    gp_9 = [gp for gp in gp_results if gp.gain_db == 9]
    if gp_9:
        K_9 = gp_9[0].K_true
        K_exp = fit.K_slope
        pct_diff = 100.0 * (K_9 - K_exp) / K_exp
        status = "PASS" if abs(pct_diff) < 5.0 else "FAIL"
        print(f"  K_true(9dB) = {K_9:.4f} DN/e-  "
              f"vs ExposureSweep K = {K_exp:.4f} DN/e-  "
              f"({pct_diff:+.1f}%)  [{status}]")

    # Check 2: n_electrons consistency across unsaturated gains
    unsat = [gp for gp in gp_results
             if gp.sat_hi < DEFAULT_SAT_THRESH and np.isfinite(gp.n_electrons)]
    if unsat:
        ne_vals = np.array([gp.n_electrons for gp in unsat])
        ne_mean = np.mean(ne_vals)
        ne_std = np.std(ne_vals)
        ne_cv = 100.0 * ne_std / ne_mean if ne_mean > 0 else float("nan")
        print(f"  n_electrons (unsaturated): mean = {ne_mean:.1f} e-, "
              f"std = {ne_std:.1f} e-, CV = {ne_cv:.1f}%  "
              f"({len(unsat)} points)")

    # Summary of read noise in electrons
    unsat_rne = [gp.read_noise_e for gp in unsat if np.isfinite(gp.read_noise_e)]
    if unsat_rne:
        rne_mean = np.mean(unsat_rne)
        rne_std = np.std(unsat_rne)
        print(f"  sigma_read_e (unsaturated): mean = {rne_mean:.2f} e-, "
              f"std = {rne_std:.2f} e-  ({len(unsat_rne)} points)")

    # GAINMODE check
    modes = set(gp.gain_mode for gp in gp_results)
    print(f"  GAINMODE values across sweep: {sorted(modes)}")
    print(f"  Exposure time: {gp_results[0].exptime_ns / 1e6:.0f} ms"
          if gp_results else "")

    # ---------------------------------------------------------------
    # Part 3: Gain Model Fitting (Session B)
    # ---------------------------------------------------------------
    print(f"\n\n{'=' * 70}")
    print("Gain Model Fit: K_true vs Gain (dB)")
    print("Model: log10(K) = a + b * gain_dB  =>  K(g) = K_0 * 10^(g/dB_per_decade)")
    print(f"{'=' * 70}")

    gm = fit_gain_model(gp_results, sat_thresh=DEFAULT_SAT_THRESH)

    print(f"\n  K_0 (K at 0 dB)      : {gm.K_0:.4f} DN/e-")
    print(f"  dB per decade (10x K): {gm.dB_per_decade:.1f} dB")
    print(f"  dB per doubling (2x K): {gm.dB_per_doubling:.1f} dB")
    print(f"  Scale factor          : {gm.scale_factor:.2f}x  "
          f"(measured dB / standard voltage dB)")
    print(f"  R-squared             : {gm.r_squared:.6f}")
    print(f"  Points used           : {gm.n_points_used}  "
          f"(gain range {gm.gain_range[0]:.0f} - {gm.gain_range[1]:.0f} dB)")
    print(f"  sigma_read_e (mean)   : {gm.sigma_read_e_mean:.2f} e-")

    # Model prediction at 9 dB vs exposure sweep
    K_pred_9 = gm.K_0 * 10.0 ** (9.0 / gm.dB_per_decade)
    K_exp = fit.K_slope
    pct = 100.0 * (K_pred_9 - K_exp) / K_exp
    print(f"\n  Model prediction at 9 dB: K = {K_pred_9:.4f} DN/e-")
    print(f"  ExposureSweep PTC K     : {K_exp:.4f} DN/e-  ({pct:+.1f}%)")

    # K_0 vs EMVA spec at 0 dB
    K_lcg = EMVA_SPEC["LCG"]["K_dn_per_e"]
    K_hcg = EMVA_SPEC["HCG"]["K_dn_per_e"]
    print(f"\n  K_0 vs EMVA spec at 0 dB analog gain:")
    print(f"    K_0 (extrapolated)  = {gm.K_0:.4f} DN/e-")
    print(f"    EMVA LCG (0 dB)     = {K_lcg:.3f} DN/e-  "
          f"(ratio K_0/LCG = {gm.K_0 / K_lcg:.2f}x)")
    print(f"    EMVA HCG (3 dB)     = {K_hcg:.3f} DN/e-  "
          f"(ratio K_0/HCG = {gm.K_0 / K_hcg:.2f}x)")
    print(f"    K_0 is between LCG and HCG: "
          f"{'YES' if K_lcg < gm.K_0 < K_hcg else 'NO'}"
          f"  (GAINMODE=1 is a specific CG configuration)")

    # dB calibration table
    print(f"\n  {'Gain (dB)':>16s}  {'Std Voltage dB':>15s}  "
          f"{'Predicted K':>12s}")
    print(f"  {'-' * 48}")
    for g_db in [0, 3, 6, 9, 12, 15]:
        std_db = g_db / gm.scale_factor
        K_pred = gm.K_0 * 10.0 ** (g_db / gm.dB_per_decade)
        print(f"  {g_db:>16.1f}  {std_db:>15.1f}  {K_pred:>12.4f}")


if __name__ == "__main__":
    main()
