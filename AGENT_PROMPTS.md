# Agent Prompts — Context Window Refresh

Use these prompts to quickly re-establish context when starting a new conversation or after context compression.

**IMPORTANT**: Copy the prompt text (inside the ``` blocks) and paste it as your first message to a new Claude Code agent.

**Current status: Phase 3 IN PROGRESS. Theory page done. Use Prompt 3 to continue enhancements.**

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
