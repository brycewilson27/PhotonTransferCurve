# Project Progress

## Status: Phase 3 In Progress — 6 Pages Built (Theory page added)

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
- [x] Page: Gain Sweep Analysis — metrics vs gain curves, optimal gain identification
- [x] Testing with both datasets (ExposureSweep & GainSweep)
- [x] Page: Frame Explorer — per-setting image browser with frame previews, pixel stats, histograms, difference image, and per-pair PTC stats

## Phase 3: Enhancements
- [x] Page: PTC Theory & Derivation — full mathematical walkthrough, exposure sweep vs gain sweep comparison, clamp effects
- [ ] Linearity analysis (mean signal vs exposure time)
- [ ] DSNU / PRNU spatial maps
- [ ] Per-pixel PTC histograms
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
