# Project Progress

## Status: Phase 3 In Progress — 6 Pages Built, GitHub repo live

---

## Phase 0: Discovery & Planning
- [x] Inventory all repository files (324 FITS, 1 MATLAB script, 2 datasets)
- [x] Review MATLAB code `CMOS_PhotonTransferCurve_v3_both.m` (195 lines)
- [x] Research PTC protocol from published sources (Janesick, EMVA 1288, Harvest Imaging)
- [x] Document data structure (ExposureSweep: 18 levels, GainSweep: 24 levels)
- [x] Identify sensor characteristics (12-bit native, 8-bit stored via CLAMP — scale=1, not 16.0588)
- [x] Write `PLAN.md` with 3-phase implementation roadmap
- [x] Create `LESSONS_LEARNED.md`, `PROGRESS.md`, `AGENT_PROMPTS.md`

## Phase 1: Core Python PTC Engine (`ptc_analysis.py`)
- [x] `load_fits_image(path)` — read FITS -> numpy float64 array
- [x] `make_master_bias(bias_paths)` — median-combine dark frames
- [x] `collect_flat_pairs(folder)` — pair bright/dark frames by `dark_` filename prefix, return local darks
- [x] `mean_var_from_pair(...)` — bias-subtract, ROI extract, compute mean & 0.5*var(A-B)
- [x] `apply_quantization_correction(variance)` — subtract 1/12 DN^2 floor
- [x] ~~`rescale_to_12bit`~~ — **REMOVED**: stored DN = native DN (clamp, not rescale). No conversion needed.
- [x] `fit_ptc(means, variances)` — linear regression with saturation masking
- [x] `derive_metrics(K, intercept, adu_max)` — conversion gain, read noise, FWC, DR
- [x] Validate Python output on ExposureSweep dataset — R^2 = 0.996, 9/16 points used
- [x] `run_ptc_analysis(...)` — high-level pipeline with local dark bias support

## Phase 2: Streamlit App (`app.py`)
- [x] Project scaffolding (requirements.txt, folder structure)
- [x] Sidebar: dataset selector, ROI controls, threshold sliders
- [x] Page: Sensor Overview — FITS metadata, dataset summary, sample image viewer
- [x] Page: PTC Analysis — interactive Plotly PTC plot with fit overlay
- [x] Page: Sensor Metrics Dashboard — key metric cards, comparison table, CSV export
- [x] Page: Gain Sweep Characterization — signal, dark level, noise, saturation, SNR vs gain
- [x] Testing with both datasets (ExposureSweep & GainSweep)
- [x] Page: Frame Explorer — per-setting image browser with frame previews, pixel stats, histograms, difference image, and per-pair PTC stats

## Phase 3: Enhancements
- [x] Page: PTC Theory & Derivation — full mathematical walkthrough, exposure sweep vs gain sweep comparison, clamp effects
- [x] GitHub repo created — brycewilson27/PhotonTransferCurve (public, main branch)
- [x] Streamlit Cloud deployment — auto-download from GitHub Releases, awaiting deploy at share.streamlit.io
- [ ] Linearity analysis (mean signal vs exposure time)
- [ ] DSNU / PRNU spatial maps
- [ ] Per-pixel PTC histograms
- [x] EMVA 1288 spec values integrated (LCG/HCG reference data from FRAMOS datasheets)
- [ ] EMVA 1288 compliance report format
- [ ] Multi-camera comparison support

---

## Session Log

### Session 1 — 2026-02-13
- Explored repository: 324 FITS files across 2 sweep datasets, 1 MATLAB script
- Reviewed MATLAB code in detail — understood full pipeline from bias through fit to reporting
- Researched PTC methodology from Harvest Imaging, PhotonsToPhotos, ESA Pyxel
- Created `PLAN.md` with complete implementation roadmap
- Created project documentation files (LESSONS_LEARNED, PROGRESS, AGENT_PROMPTS)
- **CRITICAL CORRECTION**: The 8-bit storage is a CLAMP (truncate at 255), NOT a rescale. 1 stored DN = 1 native DN. The MATLAB code's `s = 16.0588` rescaling is wrong for this sensor. Updated all docs to reflect this. Removed `rescale_to_12bit` from the pipeline.
- **Next**: Begin Phase 1 — Python PTC engine implementation

### Session 2 — 2026-02-13
- Implemented all Phase 1 functions in `ptc_analysis.py` (~540 lines)
- Installed dependencies: astropy, numpy, scipy
- **Discovery**: Each subfolder has 4 frames (2 bright + 2 dark), NOT 4 flats. The ordering of bright/dark varies across subfolders (3-9ms: dark first; 12ms+: bright first).
- Built auto-detection: classify frames by mean signal level, pair only bright frames together.
- **Discovery**: 10ms and 11ms subdirectories in ExposureSweep are empty (0 files). 16 of 18 exposure levels have data.
- **User request**: Use local dark pairs as per-subfolder bias instead of global BlackImage bias. Implemented `use_local_dark=True` mode.
- Validated on ExposureSweep_Gain90db: K=3.13 DN/e-, read noise=8.71 DN (2.78 e-), FWC=81 e- (clamp-limited), R^2=0.996
- 9 of 16 points pass saturation filter (3ms excluded for sat_lo=0.69%, 21-36ms excluded for sat_hi)
- Unicode fix: replaced special characters with ASCII for Windows console compatibility
- **Next**: Phase 2 — Streamlit app

### Session 3 — 2026-02-13
- Renamed all dark frames in both datasets with `dark_` prefix (160 files total: 64 ExposureSweep + 96 GainSweep)
- Simplified `collect_flat_pairs()` — now uses `dark_` filename prefix instead of mean-signal heuristic (faster, deterministic, no FITS reads needed for classification)
- Verified analysis output is identical after rename (K=3.13, R^2=0.996)
- `run_ptc_analysis()` supports `use_local_dark` toggle for comparing local vs global bias in the app
- Updated all project docs for next-agent handoff
- **Next**: Phase 2 — Streamlit app (use Prompt 3 in AGENT_PROMPTS.md)

### Session 4 — 2026-02-13
- Built complete Streamlit app (`app.py`, ~470 lines) with all 4 pages
- Created `requirements.txt` with all dependencies (astropy, numpy, scipy, streamlit, plotly, pandas)
- Installed streamlit, plotly, pandas
- **Sidebar controls**: dataset selector, ROI center/size, saturation threshold slider, quantization correction toggle, bias method radio (local dark vs global — disabled for GainSweep)
- **Sensor Overview page**: dataset info, sensor specs, FITS header sample, Plotly heatmap image viewer with red ROI overlay
- **PTC Analysis page**: interactive Plotly scatter (Var vs Mean) with used/excluded points, dashed fit line, 255 DN clamp vertical line, log-log toggle, K/R^2/points-used metrics
- **Sensor Metrics Dashboard**: metric cards (conversion gain, read noise, FWC, DR, SNR_max, K, intercept, R^2), FWC clamp caveat warning, per-pair summary table with color-coded Used column, CSV export button
- **Gain Sweep Analysis**: mean signal vs gain, variance vs gain, saturation fraction vs gain, full PTC plot for all gain levels, per-gain summary table. Prompts user to switch dataset if on ExposureSweep.
- Verified both datasets run through the pipeline: ExposureSweep K=3.13 R^2=0.996, GainSweep K=8.39 R^2=0.976
- App tested on localhost:9501 — HTTP 200 OK
- Used `@st.cache_data` on analysis with serializable dict output (avoids numpy array hashing issues)
- **Next**: Phase 3 (future enhancements) or user testing/refinement

### Session 5 — 2026-02-13
- Built **Frame Explorer** page (~200 lines), completing Phase 2 app
- Features: setting selector dropdown, side-by-side frame previews (Bright A, B, Dark) as Plotly heatmaps with ROI overlay
- Pixel statistics table: mean, std, min, max for full image and ROI across all frames
- ROI intensity histogram with selectable frame, 0/255 clamp boundary lines
- (A-B) difference image heatmap (ROI only, RdBu diverging colorscale centered at 0)
- Per-pair PTC stats card: mean signal, variance, sat_hi/lo, fit inclusion status
- Mini PTC plot with gold star marker highlighting the currently selected pair
- Bias comparison panel: local dark vs global master bias stats side-by-side (ExposureSweep only)
- Added `mean_var_from_pair` to app.py imports for bias comparison
- App verified: HTTP 200 OK on localhost:9505
- **Phase 2 now fully complete** — all 5 pages operational
- **Next**: Phase 3 (future enhancements) or user testing/refinement

### Session 6 — 2026-02-13
- Built **PTC Theory & Derivation** page (~200 lines), first Phase 3 enhancement
- 7 sections: Introduction, Signal Chain Model, Core PTC Derivation (Poisson -> variance scaling -> linear PTC equation), Pair-Difference Method derivation, Derived Metrics (with LaTeX equations), Exposure Sweep vs Gain Sweep comparison (with mathematical analysis of why gain sweep violates the single-line assumption), 255 DN Clamp Effect
- All equations rendered via `st.latex()` with proper LaTeX formatting
- Comparison table showing property-by-property differences between sweep types
- Expandable deeper derivation for "why can't we fit per-gain PTC"
- Added page to navigation as first item (educational reference)
- App running on localhost:9501 — now has 6 pages total
- **Next**: Remaining Phase 3 enhancements or user testing

### Session 7 — 2026-02-13
- Created `.gitignore` (excludes FITS data, __pycache__, IDE files, .streamlit/)
- Installed GitHub CLI (`gh`) via winget
- Authenticated `gh` as `brycewilson27` via browser device flow
- Configured git user identity (brycewilson27)
- Initialized git repo, created initial commit (8 files, 2676 insertions)
- Created public GitHub repo: https://github.com/brycewilson27/PhotonTransferCurve
- Pushed to `main` branch successfully
- Added Prompt 7 (Streamlit Cloud Deployment) to AGENT_PROMPTS.md
- **Note**: FITS data files (~807 MB) excluded from repo — Streamlit Cloud deployment will need demo/fallback mode
- **Next**: Deploy to Streamlit Cloud (use Prompt 7), or continue Phase 3 enhancements

### Session 8 — 2026-02-13
- **Streamlit Cloud deployment preparation**
- Created zip archive of all FITS data (324 files, 806 MB -> 594 MB compressed)
- Created GitHub Release `v1.0-data` with zip as asset: https://github.com/brycewilson27/PhotonTransferCurve/releases/tag/v1.0-data
- Added auto-download logic to `app.py`: `_ensure_data()` downloads + extracts FITS data on first cloud run
  - Detects missing `PhotonTransferCurveArtifacts/` directory
  - Downloads 594 MB zip from GitHub Releases with progress bar
  - Extracts to `PhotonTransferCurveArtifacts/`, cleans up zip
  - Uses `st.rerun()` to reload app after extraction
- Updated `.gitignore` to exclude `.zip` files
- Verified app still works locally (HTTP 200 on localhost:9510)
- Pushed to `main` branch
- **Next**: Complete Streamlit Cloud deployment at https://share.streamlit.io, update with live URL

### Session 9 — 2026-02-14
- **Redesigned Gain Sweep page**: removed all PTC terminology and reframed as "Gain Sweep Characterization"
- Renamed page from "Gain Sweep Analysis" to "Gain Sweep Characterization" in navigation
- New page sections: Signal Amplification (bright + dark overlay), Dark Level vs Gain, Temporal Noise vs Gain, Saturation vs Gain, SNR vs Gain, per-gain summary table with CSV export
- Added dark frame ROI mean/std extraction from demo_data.json frame_explorer data (demo mode) and live FITS (live mode)
- Added explainer for why noise drops at high gain (clamp artifact, not real noise reduction)
- **Updated Theory page** (Section 6): replaced "Exposure Sweep vs Gain Sweep" comparison with "Why Vary Exposure, Not Gain?" — shorter, clearer explanation of why gain sweep is not a valid PTC, with pointer to the new characterization page
- Removed the gain sweep PTC fit line, K metric, R^2, and "used in fit" concepts from the gain sweep page entirely
- **Next**: Continue Phase 3 enhancements or user testing

### Session 10 — 2026-02-14
- **Integrated EMVA 1288 spec values** from FRAMOS characterization reports for Sony IMX900-AMR-C
- Extracted data from two PDF datasheets: LCG (0 dB, Low Conversion Gain) and HCG (3 dB, High Conversion Gain)
- Added `EMVA_SPEC` dict to `ptc_analysis.py` (~65 lines) with all key parameters: K, 1/K, QE, read noise, FWC, DR, SNR, DSNU, PRNU, linearity error, black level, dark current
- **Sensor Overview page**: added full EMVA 1288 spec table (LCG vs HCG), sensor model/pixel info, expandable explainer on how specs relate to our 90 dB measurements
- **PTC Analysis page**: added side-by-side comparison table (Our PTC vs Spec LCG vs Spec HCG) after derived metrics section
- **Sensor Metrics Dashboard**: added EMVA spec comparison metric cards with delta annotations
- **Theory page**: enhanced Signal Chain Model section with explicit gain stage in signal chain diagram, added expandable "Gain architecture of the IMX900" panel explaining two-layer gain (CG mode + analog gain)
- **Key insight confirmed**: K=3.13 DN/e- at 90 dB is consistent with the spec's K=0.40 (LCG) and K=1.67 (HCG) — higher analog gain produces proportionally higher K. The 90 dB gain register mapping is proprietary (not standard voltage dB).
- Updated MEMORY.md with EMVA spec data, sensor model identification, gain architecture notes
- **Next**: Commit and push; continue Phase 3 enhancements or frame simulation work

### Session 11 — 2026-02-14 (Gain Investigation Session A)
- **Added per-point gain sweep analysis** to `ptc_analysis.py` (~120 lines)
- New `GainPointResult` dataclass: gain_db, mean_signal, pair_variance, dark_pair_variance, K_apparent, K_true, read_noise_dn, read_noise_e, n_electrons, sat_hi, sat_lo, gain_mode, exptime_ns
- New `analyze_gain_sweep_per_point()` function: computes per-gain-level K estimates with read-noise correction
- K_apparent (naive var/mean) is biased high by ~28% at 90 dB due to read noise; K_true uses quadratic correction
- **Key results (GainSweep, 24 gain levels 30-220 dB):**
  - K_true(90 dB) = 3.19 DN/e- vs ExposureSweep K = 3.13 DN/e- (+1.9%) -- **PASS** (<5% threshold)
  - n_electrons = 28.0 +/- 1.2 e- (CV=4.2%) across 16 unsaturated points -- constant illumination confirmed
  - sigma_read_e = 1.39 +/- 0.15 e- (input-referred, approximately constant across gains)
  - GAINMODE = 1 at all 24 gain levels (single CG mode throughout sweep)
  - Saturation onset at ~150 dB (5% sat_hi); gains 180+ dB are heavily clamp-distorted
- Extended CLI test with verification checks and per-gain summary table
- **Next**: Session B (gain model fitting) or Session C (app integration)

### Session 12 — 2026-02-14 (Gain Investigation Session B)
- **Added gain model fitting** to `ptc_analysis.py` (~60 lines new code)
- New `GainModelResult` dataclass: K_0, dB_per_decade, dB_per_doubling, scale_factor, r_squared, n_points_used, gain_range, sigma_read_e_mean
- New `fit_gain_model()` function: fits log10(K_true) = a + b * gain_dB via scipy.stats.linregress
- **Key results (16 unsaturated points, 30-140 dB):**
  - K_0 = 1.19 DN/e- (between EMVA LCG=0.40 and HCG=1.67 -- GAINMODE=1 is a specific CG configuration)
  - dB_per_decade = 212.5 (AstroTracker dB for 10x K increase)
  - dB_per_doubling = 64.0 (AstroTracker dB for 2x K increase)
  - scale_factor = 10.63x (AstroTracker dB / standard voltage dB of 6.02)
  - R^2 = 0.996 (excellent log-linear fit)
  - Model prediction at 90 dB: K = 3.16 DN/e- vs ExposureSweep K = 3.13 (+0.8%)
  - sigma_read_e_mean = 1.39 e- across all unsaturated gain levels
- Extended CLI test: model parameters, 90 dB prediction, EMVA K_0 comparison, dB calibration table
- **Key insight**: AstroTracker gain register uses a proprietary dB scale where 64 dB doubles K (vs 6 dB standard). The scale factor of ~10.6x means 90 AstroTracker dB ~ 8.5 standard voltage dB.
- **Next**: Session C (Streamlit app integration) or further analysis

### Session 13 -- 2026-02-14 (Gain Investigation Session C -- App Integration)
- **Integrated per-point gain analysis and gain model into Streamlit app and demo data**
- Updated `generate_demo_data.py`: calls `analyze_gain_sweep_per_point()` and `fit_gain_model()`, stores per_point_analysis (24 entries) and gain_model dict in demo_data.json
- Regenerated `demo_data.json` (38.6 MB) with new gain analysis fields
- Added `cached_gain_per_point()` and `get_gain_analysis()` to app.py for live/demo mode support
  - Live mode passes ExposureSweep read_noise_e (2.78 e-) as sigma_read_e_anchor
- **4 new sections on Gain Sweep Characterization page:**
  1. **System Gain (K) vs Analog Gain** -- K_apparent (gray), K_true (blue), model fit line (orange dashed), ExposureSweep K reference star (gold). Log-scale Y axis. Metric cards: K_0, dB/decade, dB/doubling, R^2. Expander explaining quadratic correction.
  2. **Read Noise vs Analog Gain** -- DN (growing) and e- (flat) traces with mean annotation. Expander on why input-referred noise is constant.
  3. **Electron Count Consistency Check** -- n_electrons scatter with mean +/- std band. CV metric. Key self-consistency validation.
  4. **dB Scale Calibration** -- Table mapping AstroTracker dB to standard voltage dB and predicted K, with EMVA LCG/HCG reference rows. Expander on proprietary dB scale.
- **Enhanced per-gain summary table** with 5 new columns: K_apparent, K_true, Read Noise (DN), Read Noise (e-), n_electrons. NaN values shown as "-".
- **Theory page Section 8**: "Gain Sweep Analysis Method" -- read-noise bias problem, quadratic correction derivation, exponential model fit, derived parameters, AstroTracker dB interpretation. All equations in st.latex() with \mathrm{}.
- Validated all results match Sessions A/B: K_true(90dB)=3.19 (+1.9% vs 3.13), n_e=28.0+/-1.2, sigma_read_e=1.39, R^2=0.996
- App verified: HTTP 200, no console errors, all data paths validated
- **Next**: Commit and push (triggers Streamlit Cloud redeploy)
