# Agent Prompts — Context Window Refresh

Use these prompts to quickly re-establish context when starting a new conversation or after context compression.

**IMPORTANT**: Copy the prompt text (inside the ``` blocks) and paste it as your first message to a new Claude Code agent.

**Current status: Phase 3 IN PROGRESS. Sessions A and B COMPLETE (per-point K analysis + gain model fitting). Use Prompt 8c for Session C (app integration), Prompt 7 for Streamlit Cloud deployment, or Prompt 3 for general enhancements.**

---

## Prompt 1: Full Project Context Reload

Use this when resuming work and you're not sure which phase you're in — it reads everything and picks up where you left off.

```
I'm building a Photon Transfer Curve (PTC) analysis Streamlit app for CMOS sensor characterization.

Working directory: C:\Users\Bryce Wilson\Documents\PhotonTransferCurve

**Read these files first** to understand the project:
1. `PROGRESS.md` — what's done and what's next (READ THIS FIRST)
2. `LESSONS_LEARNED.md` — technical pitfalls and decisions (especially "Deep Dive: Effects of the 255 DN Clamp")
3. `PLAN.md` — full implementation plan (3 phases)
4. `ptc_analysis.py` — the Python PTC engine (Phase 1 COMPLETE)
5. `app.py` — the Streamlit app (6 pages, Phase 3 in progress)

**Repository structure:**
- `ptc_analysis.py` — Phase 1 engine, fully validated (~555 lines)
- `app.py` — Streamlit app, 6 pages (~1100 lines)
- `PhotonTransferCurveArtifacts/ExposureSweep_Gain90db/` — 16 usable exposure levels (10ms/11ms empty), 2 bright + 2 dark frames per subfolder, global bias in root
- `PhotonTransferCurveArtifacts/GainSweep_10msExpoTime/` — 24 gain levels (30-220dB), 2 bright + 2 dark frames per subfolder, no global bias
- Dark frames prefixed with `dark_`. Bright frames are `astrotracker_pic*.fits`.
- Sensor: CMOS AstroTracker, native 12-bit ADC, stored as 8-bit via CLAMP (not rescale)
- CRITICAL: 1 stored DN = 1 native DN. Values 0-255 pass through, 256-4095 clip to 255. No scale factor.

**Key math:** Var(signal) = K * Mean(signal) + sigma^2_read
- K = DN/e- (slope), 1/K = conversion gain (e-/DN)
- Pair variance: 0.5 * var(A - B) cancels fixed-pattern noise
- Quantization correction: subtract 1/12 DN^2 (12-bit ADC step, same as stored DN)
- ADU max = 255 (clamp ceiling). FWC measured here is "through the clamp window", not true sensor FWC.
- `use_local_dark=True` uses per-subfolder dark pairs as bias; `False` uses global BlackImage bias

**MANDATORY: Documentation maintenance**
As you work, you MUST keep the project docs updated for traceability and future agents:
- `PROGRESS.md`: Check off tasks as completed. Add session log entries.
- `LESSONS_LEARNED.md`: Add new insights, bugs, or decisions discovered during work.
- `PLAN.md`: Update if the plan needs to change.
- `AGENT_PROMPTS.md`: Update if prompts need more context or corrections.
Do NOT batch these — update immediately after each milestone.

Review PROGRESS.md and continue from where the last agent left off.
```

---

## Prompt 2: Resume Phase 1 — Python Engine

**COMPLETED** — Phase 1 is done. Use Prompt 3 (Phase 3) for next steps.

---

## Prompt 2b: Resume Phase 2 — Streamlit App

**COMPLETED** — Phase 2 is fully done (all 5 pages). Use Prompt 3 (Phase 3) for next steps.

---

## Prompt 3: Phase 3 — Enhancements

**Use this prompt now. This is the next step.**

```
I'm extending a Photon Transfer Curve (PTC) analysis Streamlit app with Phase 3 enhancements.

Working directory: C:\Users\Bryce Wilson\Documents\PhotonTransferCurve

**Read these files first:**
1. `PROGRESS.md` — checklist of what's done (READ THIS FIRST)
2. `LESSONS_LEARNED.md` — technical context (especially clamp effects, pitfalls, and local vs global bias)
3. `PLAN.md` (Phase 3 section) — enhancement roadmap
4. `ptc_analysis.py` — the Python engine (COMPLETE — read it to understand the API)
5. `app.py` — the Streamlit app (6 pages — read it to understand current UI)

**What's already built:**
- `ptc_analysis.py` (~555 lines): Core PTC engine with `run_ptc_analysis()`, `collect_flat_pairs()`, `mean_var_from_pair()`, `fit_ptc()`, `derive_metrics()`, dataclasses (PairResult, PTCFitResult, FlatPairInfo)
- `app.py` (~1100 lines): Streamlit app with 6 pages:
  1. PTC Theory & Derivation — full math walkthrough, exposure sweep vs gain sweep comparison, clamp effects
  2. Sensor Overview — FITS metadata, dataset summary, sample image with ROI overlay
  3. PTC Analysis — interactive Plotly PTC plot (Var vs Mean), linear fit overlay, log-log toggle
  4. Sensor Metrics Dashboard — metric cards, per-pair table, CSV export, FWC clamp caveat
  5. Gain Sweep Analysis — metrics vs gain curves, overall sweep PTC fit
  6. Frame Explorer — per-setting image browser, frame previews, pixel stats, histogram, (A-B) diff image, mini PTC with highlighted point, bias comparison

**Engine API summary (`ptc_analysis.py`):**
- `run_ptc_analysis(sweep_folder, bias_paths, roi, sat_thresh, apply_quant_corr, use_local_dark)` -> `(PTCFitResult, list[PairResult])`
- `PTCFitResult` fields: K_slope, intercept, r_squared, conversion_gain, read_noise_dn, read_noise_e, fwc_e, dynamic_range_db, snr_max, n_points_total, n_points_used, fit_mask
- `PairResult` fields: path_a, path_b, mean_signal, variance, sat_hi, sat_lo, label
- `collect_flat_pairs(folder)` -> `list[FlatPairInfo]` (bright_a, bright_b, dark_a, dark_b, label)
- `mean_var_from_pair(path_a, path_b, master_bias, roi)` -> `PairResult`
- `make_master_bias(bias_paths)` -> `np.ndarray`
- Constants: ADU_MAX=255, DEFAULT_ROI=(932,676,200,200), DEFAULT_SAT_THRESH=0.005

**Phase 3 enhancement candidates (from PROGRESS.md):**
- [x] PTC Theory & Derivation page (full math walkthrough, exposure vs gain sweep comparison)
- [ ] Linearity analysis (mean signal vs exposure time)
- [ ] DSNU / PRNU spatial maps
- [ ] Per-pixel PTC histograms
- [ ] EMVA 1288 compliance report format
- [ ] Multi-camera comparison support

**Data structure:**
- `PhotonTransferCurveArtifacts/ExposureSweep_Gain90db/` — 16 usable exposure levels (10ms/11ms empty), 2 bright + 2 dark per subfolder
  - Subfolder names encode exposure time: "3ms", "4ms", ..., "36ms"
  - Dark frames prefixed with `dark_`. Bright frames are `astrotracker_pic*.fits` (no prefix).
  - Global bias: `BlackImage.fits`, `Black_Image2.fits` in the root
- `PhotonTransferCurveArtifacts/GainSweep_10msExpoTime/` — 24 gain levels (30-220dB), same bright/dark structure
  - Subfolder names encode gain: "30db", "40db", ..., "220db"
  - No global bias — local darks only

**Streamlit app architecture:**
- Single `app.py` with sidebar radio navigation (6 pages)
- `@st.cache_data` on analysis — returns plain dicts (not dataclasses) for hashability
- Sidebar: dataset selector, ROI controls, saturation threshold slider, quantization correction toggle, bias method radio
- All analysis in native DN (stored 8-bit = native 12-bit, clamped at 255). NO rescaling.

**Critical rules:**
- All analysis in native DN. NO rescaling. 1 stored DN = 1 native DN.
- FWC must ALWAYS display caveat: "clamp-limited at 255 DN, not true sensor FWC"
- Windows console: use ASCII only in print statements (no Unicode superscripts)
- Use only uncompressed FITS files (no `_compressed` files)
- Use pathlib for all file paths (Windows compatibility)

**Dependencies:** astropy, numpy, scipy, streamlit, plotly, pandas (in requirements.txt)

**MANDATORY: Documentation maintenance**
As you work, keep project docs updated:
- `PROGRESS.md`: Check off tasks. Add session log entry.
- `LESSONS_LEARNED.md`: Add new insights or bugs.
- `PLAN.md`: Update if plan changes.
- `AGENT_PROMPTS.md`: Update if prompts need corrections.

Continue implementing the next unchecked item in PROGRESS.md Phase 3.
```

---

## Prompt 4: Debug / Fix Issues

```
I'm working on a PTC analysis Streamlit app and need to fix an issue.

Working directory: C:\Users\Bryce Wilson\Documents\PhotonTransferCurve

**Read these files for context:**
1. `PROGRESS.md` — current status
2. `LESSONS_LEARNED.md` — known pitfalls (check these FIRST before debugging)
3. The specific file(s) with the issue

**Common issues to check:**
- DO NOT apply any 12-bit rescaling — stored DN = native DN (clamp, not rescale)
- Quantization correction is 1/12 DN^2 (small, ~0.083 DN^2)
- Pair variance halving done twice (once in function, once in caller)
- Compressed vs uncompressed FITS files accidentally mixed (only use non-`_compressed` files)
- Dark frames must be prefixed with `dark_` — `collect_flat_pairs()` depends on this
- GainSweep has no global bias — must use local darks or ExposureSweep bias as fallback
- Saturation at DN=255 means true value could be 255-4095 (information lost)
- Windows cp1252 console encoding — no Unicode superscripts/subscripts in print() calls
- 10ms and 11ms ExposureSweep subdirectories are empty (0 files)
- Streamlit `@st.cache_data` needs serializable returns — no raw numpy arrays or dataclasses

**MANDATORY: Documentation maintenance**
After fixing the issue:
- `PROGRESS.md`: Add a note about the fix in the session log.
- `LESSONS_LEARNED.md`: If this was a non-obvious bug, document it so future agents don't repeat it.
Do NOT skip this — traceability matters.

Describe the issue: [paste error or describe unexpected behavior]
```

---

## Prompt 5: Quick Task (No Full Context Needed)

```
Working directory: C:\Users\Bryce Wilson\Documents\PhotonTransferCurve

Quick context: PTC analysis app, CMOS sensor (12-bit ADC, 8-bit stored via CLAMP at 255,
1 DN stored = 1 DN native, no rescaling), two datasets (ExposureSweep and GainSweep),
Python+Streamlit app with 6 pages (PTC Theory, Sensor Overview, PTC Analysis, Sensor Metrics, Gain Sweep, Frame Explorer).

After completing the task, update `PROGRESS.md` session log with what you did.
If you learned something new, add it to `LESSONS_LEARNED.md`.

Task: [describe specific task]
```

---

## Prompt 6: Understanding the Clamp and Its Effects

```
I need to understand how the 255 DN clamp affects PTC analysis for this sensor.

Working directory: C:\Users\Bryce Wilson\Documents\PhotonTransferCurve

Read the "Deep Dive: Effects of the 255 DN Clamp on PTC Analysis" section in `LESSONS_LEARNED.md`.

Key context: The AstroTracker sensor has a 12-bit ADC (0-4095) but stores data clamped to 8-bit (0-255).
This is NOT a rescale — values above 255 are hard-clipped to 255. 1 stored DN = 1 native DN.

This affects:
- PTC validity (only Zone 1 data where mean << 255 is usable)
- FWC measurement (clamp-limited, not true sensor FWC)
- Dynamic range (limited by the 255 DN ceiling)
- Gain sweep behavior (high gains saturate faster, fewer usable points)
- Exposure sweep behavior (long exposures hit the clamp)

If this Q&A reveals new understanding, update `LESSONS_LEARNED.md` with the insight.

Question: [ask your specific question about clamp effects]
```

---

## Prompt 7: Deploy to Streamlit Cloud

```
I need to deploy a Photon Transfer Curve (PTC) analysis Streamlit app to Streamlit Cloud.

Working directory: C:\Users\Bryce Wilson\Documents\PhotonTransferCurve

**Read these files first:**
1. `PROGRESS.md` — current status (READ THIS FIRST)
2. `app.py` — the Streamlit app entry point
3. `requirements.txt` — dependencies

**GitHub repo:** https://github.com/brycewilson27/PhotonTransferCurve
- Branch: `main`
- Already pushed with: app.py, ptc_analysis.py, requirements.txt, .gitignore, and all .md docs
- FITS data files are excluded via .gitignore (too large for GitHub)

**Deployment info:**
- GitHub username: `brycewilson27`
- Git credentials: configured via Windows credential manager
- GitHub CLI (`gh`): installed at `C:\Program Files\GitHub CLI\gh.exe`, authenticated
- Entry point: `app.py`
- Dependencies: astropy, numpy, scipy, streamlit, plotly, pandas (in requirements.txt)

**IMPORTANT: Data files are NOT in the repo.**
The FITS data files (~807 MB) are excluded from GitHub. The Streamlit Cloud deployment will NOT have access to the raw FITS data. Options to handle this:
1. **Demo mode**: Add sample/precomputed data so the app works without FITS files (recommended)
2. **Cloud storage**: Upload FITS to S3/GCS and download at runtime
3. **Git LFS**: Use Git Large File Storage (has bandwidth limits on free tier)
4. **Graceful fallback**: App detects missing data and shows a message explaining this is a local-data tool

**Streamlit Cloud deployment steps:**
1. Go to https://share.streamlit.io
2. Sign in with GitHub account (brycewilson27)
3. Click "New app"
4. Select repo: `brycewilson27/PhotonTransferCurve`
5. Branch: `main`
6. Main file path: `app.py`
7. Click "Deploy"
8. Auto-deploys on every push to main

**Git workflow for updates:**
```bash
cd "C:/Users/Bryce Wilson/Documents/PhotonTransferCurve"
git add app.py ptc_analysis.py requirements.txt
git commit -m "description of changes"
git push origin main
# Streamlit Cloud auto-redeploys
```

**After deploying, update:**
- `PROGRESS.md`: Add deployment session log entry with the Streamlit Cloud URL
- `AGENT_PROMPTS.md`: Add the live URL to the project info

Task: [describe what you need — e.g., "set up demo mode for cloud deployment", "add graceful fallback when data is missing", "troubleshoot deployment error"]
```

---

## Prompt 8a: Gain Investigation Session A -- Per-Point K Analysis

**Use this to start the gain investigation. Session A computes per-gain K values from the gain sweep data.**

```
I'm extending a Photon Transfer Curve (PTC) analysis app with gain investigation analysis.

Working directory: C:\Users\Bryce Wilson\Documents\PhotonTransferCurve

**Read these files first:**
1. `PROGRESS.md` -- current status (READ THIS FIRST)
2. `ptc_analysis.py` -- core engine (~620 lines). Read the full file -- you need to understand the existing dataclasses (FlatPairInfo, PairResult, PTCFitResult), constants (ADU_MAX, DEFAULT_ROI, EMVA_SPEC), and functions (collect_flat_pairs, load_fits_image, mean_var_from_pair, make_master_bias, run_ptc_analysis, fit_ptc, derive_metrics).
3. `LESSONS_LEARNED.md` -- technical pitfalls

**Background:** We have two datasets:
- ExposureSweep_Gain90db: 16 exposure levels at fixed 90 dB gain. PTC fit gives K=3.13 DN/e-, read_noise_e=2.78 e-, R^2=0.996.
- GainSweep_10msExpoTime: 24 gain levels (30-220 dB) at fixed 10 ms exposure. 1 bright pair + 1 dark pair per gain level.

The gain sweep has variance and mean per gain level, but currently only computes an overall (meaningless) fit. We want per-point K estimates at each gain level to show how K changes with analog gain.

**Task: Session A -- Per-Point Gain Analysis**

Add the following to `ptc_analysis.py`:

1. **`GainPointResult` dataclass** with fields: gain_db (int), mean_signal (float, DN), pair_variance (float, DN^2), dark_pair_variance (float, DN^2 from dark pair-diff), K_apparent (float, variance/mean -- biased by read noise), K_true (float, read-noise-corrected), read_noise_dn (float, sqrt of dark_pair_variance), read_noise_e (float, read_noise_dn/K_true), n_electrons (float, mean/K_true), sat_hi (float), sat_lo (float).

2. **`analyze_gain_sweep_per_point()` function** that:
   - Takes sweep_folder, roi, sigma_read_e_anchor (default 2.78 from exposure sweep), apply_quant_corr (default True)
   - For each gain level, uses collect_flat_pairs() to get bright and dark pairs
   - Computes dark pair-diff read noise: 0.5 * var(dark_A_roi - dark_B_roi)
   - Computes K_apparent = (pair_variance - quant_corr) / mean_signal
   - Computes K_true by solving: variance = K*mean + K^2*sigma_read_e^2 (quadratic in K, take positive root)
   - Computes n_electrons = mean / K_true and read_noise_e = read_noise_dn / K_true
   - Reads FITS header from one bright frame per level to extract GAIN, GAINMODE, EXPTIME
   - Returns list[GainPointResult]

3. **Command-line test** at the bottom (inside `if __name__ == "__main__":`):
   - Run analyze_gain_sweep_per_point() on the GainSweep folder
   - Print a table of results: gain_dB, mean, var, K_apparent, K_true, read_noise_dn, read_noise_e, n_electrons
   - Verify: K_true at 90 dB should be within 5% of 3.13 (exposure sweep K)
   - Verify: n_electrons should be approximately constant (~28 e-) across unsaturated gains
   - Use ASCII-only print formatting (no Unicode -- Windows cp1252 console)

**Key math for K_true correction:**
The PTC equation at a single point is: var = K*mean + K^2 * sigma_read_e^2
Rearranging: K^2 * sigma_read_e^2 + K * mean - var = 0
Quadratic formula: K = (-mean + sqrt(mean^2 + 4*sigma_read_e^2*var)) / (2*sigma_read_e^2)
Take the positive root only. If discriminant is negative or K <= 0, mark as invalid.

**Existing functions to reuse:**
- collect_flat_pairs(folder) -> list[FlatPairInfo] with bright_a, bright_b, dark_a, dark_b, label
- load_fits_image(path) -> np.ndarray (float64)
- mean_var_from_pair(path_a, path_b, master_bias, roi) -> PairResult (already computes mean_signal, variance, sat_hi, sat_lo)
- make_master_bias(bias_paths) -> np.ndarray

**MANDATORY: Documentation maintenance**
After completing the implementation:
- `PROGRESS.md`: Add session log entry summarizing what was done and key numerical results
- `LESSONS_LEARNED.md`: Add insight about K_apparent vs K_true bias, and dark pair-diff vs PTC intercept read noise
- Update `AGENT_PROMPTS.md` header status line if needed
Do NOT batch these -- update immediately after verifying the code works.
```

---

## Prompt 8b: Gain Investigation Session B -- Gain Model Fitting

**Use this after Session A is complete. Fits an exponential model to K vs gain.**

```
I'm extending a Photon Transfer Curve (PTC) analysis app with gain model fitting.

Working directory: C:\Users\Bryce Wilson\Documents\PhotonTransferCurve

**Read these files first:**
1. `PROGRESS.md` -- current status (READ THIS FIRST -- Session A should be complete)
2. `ptc_analysis.py` -- core engine. Focus on the new GainPointResult dataclass and analyze_gain_sweep_per_point() function added in Session A.
3. `LESSONS_LEARNED.md` -- technical pitfalls

**Background:** Session A added per-point K_true estimates at each of 24 gain levels. K_true grows with analog gain dB. This session fits an empirical model to quantify the relationship.

**Task: Session B -- Gain Model Fitting**

Add the following to `ptc_analysis.py`:

1. **`GainModelResult` dataclass** with fields: K_0 (float, K at 0 dB), dB_per_decade (float, AstroTracker dB per 10x K increase), dB_per_doubling (float, per 2x K), scale_factor (float, ratio of AstroTracker dB to standard voltage dB -- should be ~10.6), r_squared (float), n_points_used (int), gain_range (tuple of min/max gain dB used), sigma_read_e_mean (float, mean input-referred read noise across gains).

2. **`fit_gain_model()` function** that:
   - Takes list[GainPointResult] and sat_thresh (default 0.005)
   - Filters to unsaturated points (sat_hi < sat_thresh AND K_true > 0)
   - Fits: log10(K_true) = a + b * gain_dB via scipy.stats.linregress
   - Derives: K_0 = 10^a, dB_per_decade = 1/b, dB_per_doubling = log10(2)/b, scale_factor = dB_per_doubling / 6.02
   - Computes sigma_read_e_mean as mean of read_noise_e across used points
   - Returns GainModelResult

3. **Extend the command-line test** (inside `if __name__ == "__main__":`):
   - After Session A output, call fit_gain_model() on the GainPointResult list
   - Print model parameters: K_0, dB_per_decade, dB_per_doubling, scale_factor, R^2
   - Print model prediction at 90 dB and compare to exposure sweep K=3.13
   - Print dB calibration table: a few key AstroTracker dB values -> equivalent standard voltage dB -> predicted K
   - Compare K_0 extrapolation at 0 dB to EMVA spec LCG (0.399) and HCG (1.673)

**Expected results:**
- K_true(g) = ~1.19 * 10^(g/212.5), R^2 ~ 0.996
- dB_per_decade ~ 212.5, dB_per_doubling ~ 64, scale_factor ~ 10.6
- K_0 ~ 1.19 (between EMVA LCG and HCG -- AstroTracker GAINMODE=1 is a specific CG configuration)
- sigma_read_e_mean ~ 1.0-1.5 e-

**MANDATORY: Documentation maintenance**
After completing:
- `PROGRESS.md`: Add session log entry with model fit results
- `LESSONS_LEARNED.md`: Add insight about proprietary dB scale and K_0 vs EMVA specs
- `MEMORY.md`: Add gain model parameters (K_0, dB_per_decade, scale_factor)
```

---

## Prompt 8c: Gain Investigation Session C -- App Integration

**Use this after Sessions A and B are complete. Adds gain analysis visualization to the Streamlit app.**

```
I'm extending a Photon Transfer Curve (PTC) analysis Streamlit app with gain analysis visualization.

Working directory: C:\Users\Bryce Wilson\Documents\PhotonTransferCurve

**Read these files first:**
1. `PROGRESS.md` -- current status (Sessions A and B should be complete)
2. `ptc_analysis.py` -- core engine with new GainPointResult, GainModelResult, analyze_gain_sweep_per_point(), fit_gain_model()
3. `app.py` -- Streamlit app (~1500 lines, 6 pages). Focus on page_gain_sweep() and page_ptc_derivation().
4. `generate_demo_data.py` -- demo data generator for Streamlit Cloud deployment

**Background:** Sessions A and B added per-point K computation and gain model fitting to the engine. This session integrates those results into the Streamlit app and regenerates demo_data.json.

**Verified results from Sessions A & B (use these for validation):**
- 24 gain levels total (30-220 dB), 16 unsaturated (30-140 dB, sat_hi < 0.5%)
- K_true(90 dB) = 3.19 DN/e- vs ExposureSweep K = 3.13 (+1.9%) -- PASS
- K_apparent overestimates K_true by ~28% at 90 dB due to read-noise bias
- n_electrons = 28.0 +/- 1.2 e- (CV=4.2%) across 16 unsaturated points
- Gain model: K(g) = 1.19 * 10^(g/212.5), R^2 = 0.996
- K_0 = 1.19 DN/e- (between EMVA LCG=0.40 and HCG=1.67; ratios 2.98x and 0.71x)
- dB_per_decade = 212.5, dB_per_doubling = 64.0, scale_factor = 10.63x
- sigma_read_e_mean = 1.39 +/- 0.15 e- (dark pair-diff method, constant across gains)
- GAINMODE = 1 at all 24 gain levels (single CG configuration)
- Saturation onset at ~150 dB (5% sat_hi); gains 160+ are clamp-distorted

**Critical function signatures (for live mode calls):**
- `analyze_gain_sweep_per_point(sweep_folder, roi, sigma_read_e_anchor=2.78, apply_quant_corr=True)` -> list[GainPointResult]
  - sigma_read_e_anchor MUST come from ExposureSweep PTC intercept (2.78 e- = fit.read_noise_e)
  - In app.py, pass the exposure sweep's read_noise_e as the anchor when computing live
- `fit_gain_model(gain_points, sat_thresh=0.005)` -> GainModelResult
  - sat_thresh=0.005 filters to 16 of 24 points (30-140 dB range)

**Two different read noise values -- don't confuse them:**
- PTC intercept read noise: 8.71 DN = 2.78 e- (from ExposureSweep fit y-intercept; includes all constant noise sources)
- Dark pair-diff read noise: ~4.3 DN = ~1.39 e- at 90 dB (from 0.5*var(darkA-darkB); purer readout noise floor)
- The read_noise_dn/read_noise_e in GainPointResult are the dark pair-diff values (they grow in DN with gain, stay flat at ~1.39 e-)
- The sigma_read_e_anchor=2.78 is the PTC intercept value used for the K_true quadratic correction

**Task: Session C -- App Integration**

### Step 1: Update generate_demo_data.py
- Import analyze_gain_sweep_per_point and fit_gain_model
- In the gain sweep processing, call analyze_gain_sweep_per_point() and fit_gain_model()
- Add to demo_data.json gain_sweep section:
  - "per_point_analysis": list of dicts with gain_db, mean_signal, pair_variance, K_apparent, K_true, dark_pair_variance, read_noise_dn, read_noise_e, n_electrons, sat_hi, sat_lo
  - "gain_model": dict with K_0, dB_per_decade, dB_per_doubling, scale_factor, r_squared, n_points_used, gain_range, sigma_read_e_mean
- Run generate_demo_data.py to regenerate demo_data.json
- Verify the JSON is valid and contains the new fields

### Step 2: Add 4 new sections to page_gain_sweep() in app.py

Insert AFTER the existing "Signal Amplification" section and BEFORE "Temporal Noise vs Gain":

**Section: "System Gain (K) vs Analog Gain"**
- Load per_point_analysis from demo data (or compute live)
- Plotly scatter: K_apparent (gray markers, label "K_apparent (biased)") and K_true (blue markers) vs gain_dB
- Overlay: exponential model fit line (orange dashed, from gain_model K_0 and dB_per_decade)
- Reference marker: exposure sweep K=3.13 at 90 dB (gold star, use `marker=dict(symbol="star", size=14)`)
- Only show unsaturated points (filter by sat_hi < threshold from sidebar)
- Metric cards row: K_0, dB_per_decade, dB_per_doubling, model R^2
- Expander: "How is K_true computed?" explaining the read-noise correction:
  - K_apparent = var/mean overestimates by ~28% at 90 dB
  - PTC at single point: var = K*mean + K^2*sigma_read_e^2
  - Quadratic solution: K = (-mean + sqrt(mean^2 + 4*sigma_read_e^2*var)) / (2*sigma_read_e^2)
  - Uses sigma_read_e = 2.78 e- from ExposureSweep PTC intercept as anchor

**Section: "Read Noise vs Analog Gain"**
- Plotly scatter with two Y-axes or two traces:
  - read_noise_dn vs gain_dB (increasing -- shows amplification of fixed noise floor)
  - read_noise_e vs gain_dB (flat at ~1.39 e- -- shows constant input-referred noise)
- NOTE: these are dark pair-diff read noise values, NOT the PTC intercept value (2.78 e-)
- Annotation showing mean sigma_read_e = 1.39 e- value
- Expander: "Why is read noise in electrons constant?" -- because it's an input-referred property of the readout electronics before amplification. Also note: dark pair-diff (1.39 e-) is lower than PTC intercept (2.78 e-) because the latter includes additional constant noise sources.

**Section: "Electron Count Consistency Check"**
- Plotly scatter: n_electrons vs gain_dB for unsaturated points only
- Horizontal band showing mean +/- std (use add_hrect)
- Metric: mean n_electrons = ~28.0, std = ~1.2, CV = ~4.2%
- This is the key self-consistency check: fixed illumination (10 ms) = fixed electron count regardless of gain
- Note in text: saturated points (>150 dB) show inflated n_electrons due to clamp compression

**Section: "dB Scale Calibration"**
- DataFrame table: AstroTracker dB -> equivalent standard voltage dB -> predicted K
- Compute: std_voltage_dB = astrotracker_dB / scale_factor; predicted_K = K_0 * 10^(g/dB_per_decade)
- Show rows for 0, 30, 50, 70, 90, 110, 130, 150 dB
- Add EMVA spec K values: LCG = 0.399 at 0 dB, HCG = 1.673 at 3 dB as reference rows
- Note K_0 = 1.19 falls between LCG and HCG (ratio 2.98x / 0.71x)
- Expander: "Are these real dB?" -- AstroTracker dB is a proprietary register scale (~10.6x standard voltage dB). 64 AstroTracker dB doubles K, vs 6.02 standard voltage dB. So "90 dB" is really only ~8.5 standard voltage dB of actual amplification.

### Step 3: Enhance the existing per-gain summary table
Add columns: K_apparent (DN/e-), K_true (DN/e-), Read Noise (DN), Read Noise (e-), n_electrons
Format NaN values (from saturated points) as "-" or "N/A"
Update CSV export to include these.

### Step 4: Add Theory page section 8
Add "8. Gain Sweep Analysis Method" to page_ptc_derivation() after existing section 7:
- Explain K_apparent = var/mean is biased high by read noise (~28% at 90 dB)
- Show the PTC single-point equation: var = K*mean + K^2*sigma_read_e^2
- Show the quadratic correction formula for K_true
- Explain the exponential model K(g) = K_0 * 10^(g/c), fit via log10(K) = a + b*g
- Derived parameters: K_0 = 10^a, dB_per_decade = 1/b, dB_per_doubling = log10(2)/b
- Note AstroTracker "dB" is proprietary: ~10.6x standard voltage dB, so 64 dB to double K (not 6.02)
- Mention K_0 = 1.19 falls between EMVA LCG (0.40) and HCG (1.67) because GAINMODE=1
- Use st.latex() for equations (KaTeX -- use \mathrm{} not \text{})

### Step 5: Update documentation
- PROGRESS.md: Add session log entry with what was done
- MEMORY.md: Already has gain model parameters from Session B -- verify it's current

### Step 6: Commit and push
- git add app.py ptc_analysis.py generate_demo_data.py demo_data.json PROGRESS.md
- Commit with descriptive message
- Push to main (triggers Streamlit Cloud redeploy)

### Step 7: Verification
- Launch app: `streamlit run app.py --server.port 9501`
- Switch to Gain Sweep dataset in sidebar
- Go to "Gain Sweep Characterization" page -- verify all 4 new sections render
- Verify K_true(90 dB) ~ 3.19 matches exposure sweep K=3.13 within 5%
- Verify n_electrons plot is approximately flat at ~28 e-
- Verify read_noise_e plot is approximately flat at ~1.39 e-
- Verify model fit line passes through the K_true data points
- Verify the app also works in DEMO mode (no FITS data -- uses demo_data.json)
- Go to "PTC Theory" page -- verify new Section 8 renders with equations
- Export CSV -- verify new columns present (K_apparent, K_true, read_noise_dn, read_noise_e, n_electrons)
- Check Windows console for any Unicode errors

**Important architecture notes:**
- In demo mode, per_point_analysis and gain_model come from demo_data.json
- In live mode, compute them from FITS data using the new functions
  - CRITICAL: pass sigma_read_e_anchor from the ExposureSweep PTC fit (fit.read_noise_e = 2.78)
  - If ExposureSweep hasn't been run yet, use 2.78 as default
- Use @st.cache_data patterns consistent with existing code (return plain dicts, not dataclasses)
- GainPointResult has NaN values for saturated points -- handle gracefully in plots and tables
- All st.latex() must use \mathrm{} not \text{} for KaTeX compatibility
- Console output must be ASCII-safe (cp1252)
- Streamlit ports 8501-8504 often occupied; use 9501+

**MANDATORY: Documentation maintenance**
Update PROGRESS.md, MEMORY.md after completing. Do NOT skip.
```

---

## Tips for Context Efficiency

1. **Always read PROGRESS.md first** — it tells you exactly where things stand
2. **LESSONS_LEARNED.md prevents repeated mistakes** — read it before writing analysis code. Pay special attention to the "Deep Dive: Effects of the 255 DN Clamp" section.
3. **Read `ptc_analysis.py` before modifying the engine** — understand the API, dataclasses, and `run_ptc_analysis()` parameters
4. **Read `app.py` before modifying the UI** — understand the 6-page structure, sidebar controls, and caching strategy
5. **The MATLAB file is a reference, NOT ground truth** — it has a wrong 12-bit rescaling assumption. Use it for algorithm structure but skip the `* s` / `* s^2` rescaling.
6. **Dark frames are prefixed with `dark_`** — `collect_flat_pairs()` relies on this naming convention
7. **Update docs as you go** — PROGRESS.md, LESSONS_LEARNED.md, PLAN.md, AGENT_PROMPTS.md. Future agents depend on this.
8. **FITS files are large (~3.2MB each uncompressed)** — don't read them all; sample one per subfolder
9. **Windows environment** — Python 3.10, ASCII-only console output, use pathlib for paths
10. **Streamlit caching** — `@st.cache_data` needs plain dicts with `.tolist()` for numpy arrays, not dataclasses
11. **GitHub repo** — `brycewilson27/PhotonTransferCurve`, push with `git push origin main`. FITS data excluded via .gitignore.
12. **GitHub CLI** — installed at `C:\Program Files\GitHub CLI\gh.exe`, authenticated as `brycewilson27`. Use full path or add to PATH.
13. **Streamlit Cloud** — auto-deploys from `main` branch. The app needs FITS data locally — cloud deployment requires a demo/fallback mode.
