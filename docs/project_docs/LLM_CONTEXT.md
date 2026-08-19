# docs/project_docs/ LLM Context

## Scope

Methodology and physical-bench-facing project documents — mostly non-code
artifacts (reports, spreadsheets, plots) from the design/dimensioning side of
the project, plus one code-adjacent methodology doc that `optimization/`
depends on.

## Contents

| File | Role |
|---|---|
| `CALIBRATION_METHODOLOGY.md` | Model calibration methodology for the bench: which `CorrelationCoefficients` fields are identifiable from which observables, fixed-vs-calibrated parameter split, Bayesian calibration workflow (`calibrate_ls`/`calibrate_map`/`calibrate_mcmc`), recommended test matrix. Referenced directly by `optimization/LLM_CONTEXT.md` and is the authoritative source for calibration strategy — code lives in `optimization/calibrate.py`. |
| `Rapport de dimensionnement - Echangeur Hélium - août2025.docx` | French-language sizing report for the Helium heat exchanger (Aug. 2025). Not machine-read as part of this pass; a project deliverable, not a context file. |
| `Rapport de dimensionnement [notes modèle quasi-1D].docx` | French-language notes accompanying the quasi-1D model sizing report. Not machine-read. |
| `22aôut2025-concepts_de_brûleursHX.pptx` | Burner-HX concept slide deck (2025-08-22). Not machine-read. |
| `Bilan masses HX brûleur.xlsx` | Burner-HX mass budget spreadsheet. Not machine-read. |
| `Guide de conception brûleur-HX - 03sept2025.docx` | French-language burner-HX design guide (2025-09-03). Not machine-read. |
| `JetA_temperature_vs_OF_1bar.png`, `JetA_cp_vs_OF_1bar.png`, `JetA_k_vs_OF_1bar.png`, `JetA_rho_vs_OF_1bar.png`, and their `_zoom` variants | Plots of Jet-A combustion-product property (T, cp, k, rho) vs. O/F ratio at 1 bar. Reference plots, not read as text. |

## Notes

Only `CALIBRATION_METHODOLOGY.md` is a text-context file in the sense the
other `LLM_CONTEXT.md` files in this repo describe; the rest of this folder
is project deliverable material (Office documents, a spreadsheet, PNG plots)
that a future LLM session would only need to open on explicit request (e.g.
"what does the mass budget spreadsheet say"), not as background reading.

## TODOs (only evidenced items)

None found stated in `CALIBRATION_METHODOLOGY.md` itself beyond its own
recommended test matrix (which is a bench test plan, not a code TODO).

## Change history

No meaningful git history. `CALIBRATION_METHODOLOGY.md` is undated
internally but its content (parameter defaults, e.g. `ali_c_hi=0.325`,
`salimpour_a=0.317`) matches the values still in `input_data.py`'s
`CorrelationCoefficients` as of this writing.
