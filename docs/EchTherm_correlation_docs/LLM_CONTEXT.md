# docs/EchTherm_correlation_docs/ LLM Context

## Scope

Extracted literature/reference-screenshot material backing the shell-and-tube
("EchTherm-style") correlation set used by `main_solve_shellntube.py` and
`shellTubeProp` — Bell-Delaware shell-side heat transfer/pressure drop, and
Ravigururajan-Bergles internal corrugated/ribbed-tube correlations
(`inside_tube_choice="grooved"`). Generated 2026-07-10 from PNG documentation
screenshots plus one added source PDF.

## Contents

| File | Role |
|---|---|
| `EchTherm_shell_and_tube_modeling_extraction.md` | Engineering extraction (not verbatim) of the EchTherm shell-and-tube documentation screenshots plus the Ravigururajan-Bergles paper. Six numbered sections: (1) shell-side heat transfer (Bell-Delaware `alpha_real = alpha_ideal * Jc*Jf*Jb*Js*Jr`, Colburn-factor ideal coefficient), (2) shell-side pressure drop, (3) internal corrugated/ribbed tube correlations, (4) references extracted from the screenshots, (5) immediate comparison targets for Combustor-HX, (6) screenshot coverage audit. |
| `ravigururajan1996_extracted.md` | Raw OCR/text extraction of Ravigururajan & Bergles (1996), "Development and Verification of General Correlations for Pressure Drop and Heat Transfer in Single-Phase Turbulent Flow in Enhanced Tubes" — used to complete the friction/heat-transfer correlation that was partly cut off in one of the source screenshots. Validity ranges per the paper's own abstract: `e/d` 0.01–0.2, `p/d` 0.1–7.0, `alpha/90` 0.3–1.0, Re 5000–250,000, Pr 0.66–37.6. |

## Relationship to code

This folder is reference-only (no code). The correlations it documents are
implemented in `1Dmodel/physics/heat_transfer_correlations.py` /
`friction_correlations.py` (`dispatch_nu_tube_straight`, corrugated-tube
Vicente/Cruz-style path) and `1Dmodel/physics/bell_delaware.py`. See
`docs/context/PHYSICS_CONTEXT.md` for current wiring/calibration-knob status
of the grooved-tube path (`tube_grooved_Nu_factor`, `tube_grooved_f_factor`).

## TODOs (only evidenced items)

None found stated within these files beyond the extraction's own "immediate
comparison targets" section (a reference-vs-model comparison list, not a code
TODO).

## Change history

Both files carry an explicit generation date of 2026-07-10 in their own text
(no other dated revisions evidenced).
