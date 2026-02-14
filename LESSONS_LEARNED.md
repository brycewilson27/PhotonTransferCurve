# Lessons Learned

Ongoing log of technical insights, pitfalls, and decisions made during development.

---

## Data & Sensor Quirks

### 8-bit Clipping (NOT Rescaling) — CRITICAL
- The AstroTracker CMOS sensor is natively 12-bit (0–4095 DN) but stores FITS files as 8-bit (0–255 DN).
- **The 8-bit conversion is a CLAMP, not a rescale.** Values 0–255 pass through as-is. Values 256–4095 are clipped to 255.
- This means **1 DN in the stored file = 1 DN in the native 12-bit ADC**. There is no scale factor.
- The MATLAB code's `s = 4095/255 = 16.0588` rescaling is **incorrect** for this sensor. It artificially inflates means by 16x and variances by 258x.
- The effective usable dynamic range is only the bottom 255/4095 = **6.2%** of the sensor's native 12-bit range.
- Saturation in the stored data (DN=255) means the true 12-bit value could be anywhere from 255 to 4095 — information is lost.
- **For our Python implementation: work entirely in the 8-bit domain as stored. No rescaling needed.**
- ADU max for saturation detection = 255 (the stored ceiling).

### Quantization Noise Correction
- The sensor ADC has 12-bit resolution, so the native quantization step is 1 DN_12.
- Since the clamp preserves the bottom 8 bits unchanged, the quantization noise is that of the **12-bit ADC**, not an 8-bit requantization.
- Quantization variance = `1/12 DN²` (one 12-bit LSB step, and 1 stored DN = 1 native DN).
- This is a small correction (~0.083 DN²) and may be negligible compared to read noise, but should still be applied for rigor.
- The MATLAB code's `Vq8 = 1/12` happens to be numerically correct despite the wrong reasoning — the stored DN *is* the native DN.

### Pair-wise Variance Method
- Using `0.5 * var(A - B)` instead of direct frame variance eliminates fixed-pattern noise (FPN/PRNU).
- Critical: the MATLAB helper `mean_var_from_pair` already halves the variance — doubling the halving would underestimate noise by 2x.
- Always verify the convention of your variance function before applying the 0.5 factor.

### Saturation Masking
- The MATLAB code uses a 0.5% threshold — if more than 0.5% of ROI pixels are at `aduMax` or 0, the pair is excluded from the fit.
- Both high-end (saturated) and low-end (clipped to zero) saturation are tracked separately.
- Including saturated points in the fit will pull the slope down and underestimate conversion gain.

---

## Deep Dive: Effects of the 255 DN Clamp on PTC Analysis

This section documents how the clamp (not rescale) conversion from 12-bit to 8-bit affects every aspect of our photon transfer curve analysis. This is the single most important thing to understand about this dataset.

### What the Clamp Does Physically

```
Native 12-bit ADC output:  0  1  2  ... 254  255  256  257  ... 4095
                           |  |  |       |    |    |    |         |
Stored 8-bit value:        0  1  2  ... 254  255  255  255  ... 255
```

- The camera's ADC digitizes to 12 bits (0–4095 DN).
- Before writing to the FITS file, the firmware clamps: `stored = min(native, 255)`.
- This is a **lossy, nonlinear, irreversible** operation — once a pixel hits 255, we cannot know its true value.
- This is NOT the same as bit-shifting (`native >> 4`), which would divide by 16 and preserve the full range in coarser steps.

### Why This Matters for the PTC

The PTC assumes a **linear** relationship between signal and variance over the unsaturated range. The clamp introduces three distinct zones in our data:

```
Zone 1: CLEAN DATA (mean << 255)
  - No pixels in the ROI are hitting the clamp
  - Stored values = true values
  - PTC relationship holds perfectly: Var = K * Mean + sigma²_read
  - This is where our linear fit is valid

Zone 2: PARTIALLY CLAMPED (mean approaching 255)
  - Some pixels in the ROI are at 255 (clamped), others are not
  - The mean is pulled DOWN because clamped pixels can't exceed 255
  - The variance is SUPPRESSED because the upper tail of the distribution is crushed to 255
  - Both mean and variance are biased low → the PTC point falls below the true line
  - This is the "rollover" region where the PTC curves downward

Zone 3: FULLY CLAMPED (mean ≈ 255)
  - Most/all pixels are at 255
  - Mean ≈ 255, Variance ≈ 0
  - No useful information — all signal structure is destroyed
```

### Effect on Each Derived Parameter

#### Conversion Gain (K and 1/K)
- **If we fit only Zone 1 data**: K is accurate. The slope of Var vs Mean in the unsaturated region correctly gives the system gain in DN/e⁻ because the stored DN = native DN.
- **If Zone 2 points leak into the fit**: K will be underestimated (slope biased low) because partially-clamped points have suppressed variance. This is why the saturation mask is critical.
- **Practical impact**: With a 0.5% saturation threshold, we should catch most Zone 2 contamination.

#### Read Noise (sigma_read)
- **Unaffected by the clamp** — read noise is measured from the y-intercept of the PTC fit, which corresponds to the zero-signal regime (dark/bias frames, well below 255).
- Read noise is also measured directly from the variance of bias frames, which are typically in the 0–30 DN range, far from the clamp.
- **Caveat**: If the analog gain is set so high that even dark current pushes pixels near 255, then read noise measurement would be compromised. Check that bias frame means are well below 255.

#### Full Well Capacity (FWC)
- **This is the most affected parameter.**
- The "FWC" we measure is actually the **clamp-limited capacity**: how many electrons it takes to reach 255 DN at the current gain setting.
- Formula: `FWC_measured = 255 / K` (in e⁻), where K is DN/e⁻.
- The true sensor FWC (how many electrons the pixel can physically hold before blooming) requires 12-bit data to measure — it may be 16x higher or more.
- **What this means**: Our "FWC" is really "usable FWC within the 8-bit output window." It's still a valid and useful number for characterizing the system *as operated*, but it should not be compared to datasheet FWC specs that assume full-range readout.
- **In the Streamlit app**: Always display a caveat that FWC is "as seen through the 255 DN clamp."

#### Dynamic Range
- `DR = FWC / sigma_read` (in e⁻)
- Since FWC is clamp-limited, the measured DR is also clamp-limited.
- True sensor DR = true FWC / sigma_read, which would be higher.
- **Our measured DR represents the actual usable DR of the system as configured** — arguably the more operationally relevant number.

#### Signal-to-Noise Ratio (SNR)
- SNR = Mean / sqrt(Var) = Mean / sqrt(K * Mean + sigma²_read)
- In Zone 1, this is accurate.
- Maximum SNR occurs at the clamp ceiling (255 DN), not at the true sensor FWC.
- `SNR_max ≈ sqrt(255 / K)` for the shot-noise-limited case.

### Effect on the Gain Sweep Dataset

The GainSweep dataset varies analog gain from 30 dB to 220 dB at fixed 10 ms exposure. The clamp has **different severity** at different gains:

```
Low gain (30 dB):
  - K is small (few DN per electron)
  - Signal in DN may be low → most/all pixels in Zone 1
  - PTC analysis is clean, but FWC is clamp-limited

High gain (220 dB):
  - K is large (many DN per electron)
  - Even moderate illumination pushes pixels to 255 quickly
  - More pairs will be in Zone 2/3 → more points excluded by saturation mask
  - Fewer usable data points for the linear fit
  - FWC (in e⁻) is very small because 255 DN is reached with fewer photons
```

**Key question for Gain Sweep**: At very high gains, we may have so few unsaturated pairs that the linear fit is unreliable. The Gain Sweep page reports points used / total and the overall fit R^2. Future: add per-gain-level warnings when usable points < 3.

### Effect on the Exposure Sweep Dataset

The ExposureSweep varies exposure from 3 ms to 36 ms at fixed 90 dB gain:

```
Short exposures (3 ms):
  - Low signal → Zone 1 (clean)
  - May be near the noise floor → shot noise competes with read noise

Long exposures (36 ms):
  - High signal → likely Zone 2 or Zone 3
  - These pairs get excluded by saturation masking
  - Gives us the "rollover" visible in the PTC plot
```

**Expected PTC shape**: Linear rise from ~3 ms, then rollover as exposures push pixels toward 255. The linear fit uses only the rising portion.

### How the Clamp Differs from a Rescale — Why It Matters

| Property | Clamp (our sensor) | Rescale (what MATLAB assumes) |
|----------|-------------------|-------------------------------|
| Mapping | `stored = min(native, 255)` | `stored = native / 16` (approx) |
| DN value meaning | 1 stored DN = 1 native DN | 1 stored DN ≈ 16 native DN |
| Dynamic range coverage | Bottom 6.2% of sensor range | Full sensor range at coarser resolution |
| Information loss | Above 255 DN native — total loss | Uniform — 4 bits of precision lost everywhere |
| Noise floor | Native ADC noise (1 LSB_12) | Requantization noise (~16x larger) |
| Saturation behavior | Sharp clamp at 255 (hard clip) | Gradual approach to 255 (maps to 4095) |
| PTC fit validity | Valid below 255, meaningless above | Valid across full range (with requant correction) |
| Scale factor for analysis | **s = 1** (no conversion) | s = 16.0588 (must convert) |
| Quantization correction | 1/12 DN² (12-bit ADC step) | s²/12 ≈ 21.5 DN² (8-bit requantization) |

### Practical Recommendations for Our Implementation

1. **Aggressive saturation masking**: DONE. The sidebar slider lets users tune the threshold (default 0.5%). PTC plot and Frame Explorer both visually distinguish used vs excluded points.

2. **Visual clamp indicator**: DONE. PTC plots include a vertical dashed line at Mean=255 DN. Frame Explorer histogram shows 0 and 255 boundary lines.

3. **Per-gain usable range**: DONE (partial). The Gain Sweep page shows saturation fraction vs gain with threshold line. The per-gain summary table shows sat_hi/sat_lo percentages.

4. **FWC caveat everywhere**: DONE. The Sensor Metrics Dashboard always shows a warning that FWC is clamp-limited at 255 DN. Future: consider also computing `FWC_true_estimate = 4095 / K` as a rough upper bound.

5. **Bias frame sanity check**: Frame Explorer pixel statistics table shows dark frame means, allowing manual verification. Future: could add an automated warning if dark frame ROI mean > 200 DN.

6. **Log-log PTC toggle**: DONE. PTC Analysis page has a checkbox to switch both axes to log scale.

---

## FITS File Organization

### Naming Convention
- Files: `astrotracker_pic2026-02-12--HH-MM-SS.fits`
- Each exposure/gain setting has 4 uncompressed + 4 compressed copies.

### Frame Structure per Subfolder
- Each subfolder has 4 uncompressed FITS files: **2 bright (illuminated) + 2 dark**.
- Dark frames have been renamed with a `dark_` prefix (e.g. `dark_astrotracker_pic2026-02-12--00-33-43.fits`). 160 files renamed total across both datasets.
- The ordering of bright vs dark frames varied: at lower exposures (3-9ms) dark came first; at higher exposures (12ms+) bright came first. The `dark_` prefix now makes this unambiguous.
- `collect_flat_pairs()` uses the `dark_` filename prefix to classify (no FITS reads needed).
- The 2 dark frames per subfolder serve as **local bias** for that exposure/gain level.
- Using local dark pairs as bias is preferred over global bias because it accounts for exposure-dependent dark current.

### Local vs Global Bias Comparison
- `run_ptc_analysis(use_local_dark=True)` uses dark pair from each subfolder as bias.
- `run_ptc_analysis(use_local_dark=False)` uses global `BlackImage.fits` / `Black_Image2.fits`.
- Local dark bias gives lower read noise estimate (8.71 DN vs 10.19 DN) because it better matches the dark current at each exposure time.
- K_slope is robust to bias method (~3.13 either way) because pair-difference cancels bias.
- The Streamlit app offers a bias method toggle in the sidebar (local dark vs global master bias) and the Frame Explorer page includes a per-pair bias comparison table.

### Missing Data
- `ExposureSweep_Gain90db/10ms/` and `ExposureSweep_Gain90db/11ms/` directories exist but are empty (0 files).
- This means 16 of 18 expected exposure levels have usable data.

### Black/Bias Frames
- Global bias located in `ExposureSweep_Gain90db/` root: `BlackImage.fits`, `Black_Image2.fits`
- Image shape: 1552 x 2064 pixels, uint8
- Bias frame center ROI (200x200) means: ~12.3 DN (well below 255, safe)
- The GainSweep dataset does NOT have its own bias frames — uses local dark pairs from each gain subfolder.
- Global bias used as fallback when local darks are unavailable.

---

## Algorithm & Math Notes

### PTC Linear Model
```
Var(signal) = K * Mean(signal) + sigma²_read
```
- Slope K = DN/e⁻ (system gain, EMVA convention)
- 1/K = e⁻/DN (conversion gain, Janesick convention)
- y-intercept = read noise variance in DN²
- Read noise in e⁻ = sqrt(intercept) / K

### Full Well Capacity
- FWC = ADU_max / K (in e⁻), but practically use the last unsaturated mean value
- **Because of the 255 DN clamp, the measured FWC is the capacity up to the clamp ceiling, NOT the true sensor FWC.**
- The true sensor FWC would require unclamped 12-bit data (0–4095) to measure.
- What we measure is effectively "FWC as seen through the 8-bit window" — still useful for characterizing the system as operated.

### Dynamic Range
- DR = FWC / read_noise_e (in e⁻), often reported in dB: `20 * log10(DR)`

---

## Development Decisions

### Why Python + Streamlit (not MATLAB)
- MATLAB requires expensive licenses; Python is free and deployable
- Streamlit provides instant interactive UI with minimal code
- `astropy.io.fits` is the standard for FITS in Python
- Plotly gives interactive zoom/hover that MATLAB figures can't easily serve on the web

### ROI Approach
- MATLAB uses interactive `choose_roi` — not suitable for a web app
- Streamlit sidebar uses numeric inputs (center x/y, width, height) with red rectangle overlay on image previews
- Default: center 200x200 px at (932, 676) to match MATLAB baseline

---

## Validation Results (ExposureSweep_Gain90db, Session 2)

### PTC Fit (local dark bias, center 200x200 ROI)
- K (slope) = 3.13 DN/e- (EMVA convention)
- Conversion gain = 0.32 e-/DN (Janesick convention)
- Read noise = 8.71 DN = 2.78 e- rms
- FWC (clamp-limited) = 81 e- at 255 DN ceiling
- Dynamic range = 29.3 dB
- R^2 = 0.996 (excellent linear fit)
- 9 of 16 points used (3ms excluded for floor saturation, 21-36ms excluded for ceiling saturation)

### Physical Plausibility Check
- At 90 dB analog gain, K ~3 DN/e- is expected (very high gain = many DN per electron)
- Read noise of 2.78 e- is reasonable for a CMOS sensor at high gain
- FWC of 81 e- is very small but expected: 255 DN / 3.13 DN/e- = 81 e-
- The low FWC is because the 8-bit clamp limits usable range, not because the sensor well is small
- The PTC shows textbook behavior: linear rise (4-18ms), rollover (21-24ms), full clamp (30-36ms)

---

## Pitfalls to Avoid

1. **There is only one domain** — the stored 8-bit values ARE the native 12-bit values (clamped at 255). Do NOT apply any scale factor. The MATLAB code's 12-bit rescaling is wrong for this sensor.
2. **Don't use compressed FITS** for analysis — stick to uncompressed files to avoid compression artifacts in noise measurements
3. **Don't forget bias subtraction** — raw mean includes dark current offset; subtracting master bias is mandatory before computing signal mean
4. **Don't assume gain labels are linear** — the dB gain labels (30–220 dB) may not correspond linearly to actual analog gain; verify with the PTC results
5. **Don't include edge pixels in ROI** — vignetting and optical falloff at edges will introduce spatial variance that's not sensor noise
6. **Dark frames are prefixed with `dark_`** — all 160 dark frames across both datasets have been renamed. `collect_flat_pairs()` relies on this prefix. If new data is added, dark frames must follow this convention.
7. **Windows console can't print Unicode** — use ASCII only in print statements (e.g. `e-` not `e⁻`, `DN2` not `DN²`). The cp1252 encoding will crash on special characters.
8. **Streamlit `@st.cache_data` requires serializable return types** — dataclasses with `np.ndarray` fields (like `PTCFitResult.fit_mask`) can't be hashed. Convert to plain dicts with `.tolist()` for arrays before returning from cached functions.
9. **GainSweep has only 1 bright pair per gain level** — can't fit a per-gain PTC. Instead, run the overall sweep analysis (which fits across all gain levels combined) and display per-gain-level statistics (mean, variance, saturation).
10. **Streamlit port conflicts on Windows** — multiple Streamlit instances may occupy ports 8501-8504. Use `netstat -ano | findstr ":850"` to find conflicts, or pick a higher port like 9501.
11. **Exposure sweep PTC fit is rigorous; gain sweep fit is an approximation** — In an exposure sweep, K and sigma_read are constant across all data points, so the linear model `Var = K*Mean + b` is correct. In a gain sweep, both K and sigma_read change at every point (`Var_i = K_i*Mean_i + b_i`), violating the single-line assumption. The gain sweep overall fit gives an approximate average K but is not physically rigorous. See the "PTC Theory & Derivation" page in the app for the full mathematical comparison.
12. **Streamlit `st.latex()` renders KaTeX** — use standard LaTeX math notation. Inline LaTeX in `st.markdown()` uses `$...$` syntax. For display equations use `st.latex()` which renders centered. Avoid `\text{}` with special characters; use `\mathrm{}` instead.
