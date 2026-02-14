# Photon Transfer Curve — Streamlit App Plan

## 1. Repository Inventory

### Files
- **1 MATLAB script**: `PhotonTransferCurveArtifacts/CMOS_PhotonTransferCurve_v3_both.m`
- **324 FITS files** (~807 MB total) across two datasets:
  - **ExposureSweep_Gain90db/** — 18 exposure times (3–36 ms) at fixed 90 dB gain
  - **GainSweep_10msExpoTime/** — 24 gain levels (30–220 dB) at fixed 10 ms exposure
- Each subdirectory: 4 uncompressed + 4 compressed `.fits` files
- Black/bias frames in `ExposureSweep_Gain90db/` (BlackImage, Black_Image2)
- **Sensor**: CMOS AstroTracker star tracker, native 12-bit ADC, stored as 8-bit via **clamp** (not rescale — values 0–255 pass through, 256–4095 clip to 255)
- **Important**: 1 stored DN = 1 native DN. No scale factor. Effective dynamic range is bottom 6.2% of sensor.

### MATLAB Code Analysis (`CMOS_PhotonTransferCurve_v3_both.m`)
The script performs:
1. **Master bias** from dark frames (`make_master_bias`)
2. **ROI selection** (200×200 px) on a sample flat
3. **Pair-wise variance** via `mean_var_from_pair` — computes `0.5 * var(A - B)` to cancel fixed-pattern noise
4. **Quantization correction** — subtracts 1/12 DN² floor
5. ~~**8-bit → 12-bit rescaling** — mean × s, variance × s²~~ **NOTE: This step is INCORRECT for our sensor. The 8-bit data is clamped, not rescaled. Our Python port skips this entirely.**
6. **Saturation masking** — excludes pairs where >0.5% of ROI pixels are clipped
7. **Linear fit** via `fit_ptc` → slope = K (DN/e⁻), intercept = σ²_read
8. **Derived metrics**: conversion gain (e⁻/DN), read noise (DN & e⁻), dark noise, FWC
9. ~~**Dual reporting** in 12-bit and 8-bit-equivalent domains~~ **Single domain (8-bit = native) in our port**
10. **PTC plots** — Variance vs Mean with fit line, clipped-point markers

---

## 2. PTC Protocol (from literature + MATLAB code)

### Theory
The PTC relates **mean signal** (DN) to **temporal noise variance** (DN²):

```
Var(signal) = K · Mean(signal) + σ²_read
```

Where:
- **K** = system gain in DN/e⁻ (slope of the linear region)
- **1/K** = conversion gain in e⁻/DN
- **σ²_read** = read noise variance (y-intercept)
- **FWC** = full well capacity in e⁻ (signal at saturation / K)

### Noise Regimes (log-log PTC)
| Region | Slope | Dominant Noise |
|--------|-------|----------------|
| Floor  | 0     | Read noise (thermal + electronics) |
| Linear | 0.5   | Photon shot noise (√N) |
| Rolloff| 1.0   | Fixed-pattern noise / saturation |

### Data Acquisition Protocol
1. Capture **dark/bias** frames (no light, minimum exposure)
2. For each illumination level, capture **pairs** of identically-exposed flat-field frames
3. Compute per-pair: mean of `(A+B)/2`, variance from `0.5 * var(A-B)` (cancels FPN)
4. Repeat across full dynamic range (vary exposure or gain)

---

## 3. Implementation Plan

### Phase 1: Core Python PTC Engine — COMPLETE

`ptc_analysis.py` (~555 lines) implements the full pipeline:

| Step | Function | Status | Notes |
|------|----------|--------|-------|
| 1 | `load_fits_image(path)` | Done | Returns float64 array |
| 2 | `make_master_bias(bias_paths)` | Done | Pixel-wise median combine |
| 3 | `collect_flat_pairs(folder)` | Done | Uses `dark_` filename prefix to classify bright/dark pairs |
| 4 | `mean_var_from_pair(pathA, pathB, master_bias, roi, adu_max)` | Done | Returns `PairResult` dataclass with mean, variance, sat fractions |
| 5 | `apply_quantization_correction(variance)` | Done | Subtracts 1/12 DN^2, floors at 0 |
| 6 | `fit_ptc(means, variances, sat_hi, sat_lo, sat_thresh)` | Done | OLS with saturation masking |
| 7 | `derive_metrics(K_slope, intercept, ...)` | Done | Returns `PTCFitResult` with all sensor metrics |
| 8 | `run_ptc_analysis(sweep_folder, bias_paths, ...)` | Done | High-level pipeline with `use_local_dark` toggle |

**Key deviations from original plan:**
- Added `FlatPairInfo` dataclass to carry bright + dark pairs per subfolder
- Dark frames renamed with `dark_` prefix (160 files) for deterministic classification
- Added `use_local_dark` parameter to toggle between per-subfolder dark bias and global bias
- `fit_ptc` and `derive_metrics` were separated for reuse in the Streamlit app

### Phase 2: Streamlit App (`app.py`) — FULLY COMPLETE
Built as a single `app.py` (~1100 lines) with 6 pages:

#### Sidebar
- Dataset selector (ExposureSweep vs GainSweep)
- ROI controls (center x/y, width, height) with preview
- Saturation threshold slider
- Quantization correction toggle
- **Bias method toggle**: local dark pairs (per subfolder) vs global master bias (BlackImage.fits)
- Note: all analysis in native DN (stored 8-bit = native 12-bit low byte, clamped at 255)

#### Page: PTC Theory & Derivation
Educational reference page with full mathematical walkthrough:
- Introduction to PTC concept (sensor fingerprint)
- Signal chain model (photon -> electron -> DN, gain K, noise sources table)
- Core PTC derivation from Poisson statistics through variance scaling to the linear PTC equation
- Pair-difference method derivation (canceling FPN via 0.5*var(A-B))
- Derived metrics with LaTeX equations (K, g, read noise, FWC, DR, SNR_max, quantization correction)
- **Exposure Sweep vs Gain Sweep comparison**: detailed side-by-side analysis of mathematical validity, what varies, what's constant, why gain sweep violates single-line assumption, comparison table
- 255 DN clamp effect summary (three zones, FWC caveat)
- References (Janesick, EMVA 1288, Bohndiek, Harvest Imaging)

#### Page: Sensor Overview
- Two-column layout: dataset info (path, levels, frame counts, labels) and sensor specs (ADC, format, image size)
- FITS header sample from first bright frame (up to 15 key cards)
- Sample image viewer (Plotly heatmap of middle exposure/gain level) with red ROI rectangle overlay

#### Page: PTC Analysis
- Interactive PTC plot (Variance vs Mean) — Plotly for zoom/hover
- Linear fit overlay with equation displayed in legend
- Clipped (red x) vs used (blue dots) points clearly marked
- Vertical dashed line at 255 DN clamp ceiling
- Log-log toggle checkbox to switch axes to log scale
- K, R^2, and points-used metric cards

#### Page: Sensor Metrics Dashboard
- Key metrics cards row 1: Conversion Gain (e-/DN), Read Noise (e- rms), Dynamic Range (dB), SNR_max
- FWC metric with mandatory clamp-limited caveat warning
- Key metrics cards row 2: K (slope), Intercept, Read Noise (DN), R-squared
- Per-pair summary table with color-coded "Used in Fit" column (green/red)
- CSV export button

#### Page: Gain Sweep Analysis
- Overall sweep PTC fit across all gain levels (K, R^2, points used)
- Mean signal vs analog gain plot (with 255 DN clamp line)
- Variance vs analog gain plot (color-coded by used/excluded)
- Saturation fraction vs gain (hi and lo, with threshold line)
- Full PTC scatter plot for all gain levels with fit overlay
- Per-gain summary table
- Note: GainSweep has only 1 pair per gain level, so per-gain PTC fitting is not possible — uses overall sweep fit instead

#### Page: Frame Explorer
Per-subfolder image browser with per-pair PTC statistics combined into one view:
- **Setting selector**: dropdown to pick any exposure time or gain level
- **Frame previews**: side-by-side bright A, bright B, and local dark frames as Plotly heatmaps with red ROI overlay
- **Pixel statistics table**: mean, std, min, max for each frame (full image and ROI columns)
- **ROI intensity histogram**: selectable frame (Bright A/B/Dark), 256-bin histogram with vertical clamp boundary lines at 0 and 255
- **Difference image**: (A - B) ROI heatmap with RdBu diverging colorscale (centered at 0), showing temporal noise spatially with mean/std in title
- **PTC pair stats card**: metric cards for mean signal, variance, sat_hi%, sat_lo%, plus green/red fit inclusion indicator
- **Mini PTC plot**: all PTC points shown (used=lightblue, excluded=lightcoral) with fit line and gold star marker highlighting the currently selected pair
- **Bias comparison**: if ExposureSweep with local darks and global bias both available, shows side-by-side table of per-pair stats under each bias method

### Phase 3: Enhancements — IN PROGRESS
- [x] PTC Theory & Derivation page (full math walkthrough, exposure vs gain sweep comparison)
- [ ] Linearity analysis (mean signal vs exposure time)
- [ ] DSNU / PRNU maps
- [ ] Per-pixel PTC histograms
- [ ] EMVA 1288 compliance report
- [ ] Multi-camera comparison
