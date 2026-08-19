# docs/ LLM Context

## Scope

`docs/` holds everything that is not code: canonical per-topic context files,
methodology notes, the active solver-architecture design docs, literature
source material backing the liquid/boiling correlations, and generated
validation-run output. This file is a map of the subfolders — read the
subfolder's own `LLM_CONTEXT.md` for real content.

## Contents

| Path | Role |
|---|---|
| `README.md` | Short top-level doc index (predates most of the subfolders below). |
| `TECHNICAL_REFERENCE.md` | Master technical reference: physics, numerics, materials, assumptions, validation status. Large (~61 KB), organized by numbered section (see `CLAUDE.md`'s canonical list). |
| `USER_GUIDE.md` | User-facing how-to-run guide. |
| `context/` | Canonical per-topic context files explicitly listed in root `CLAUDE.md` — solver architecture, physics, transient status, and one comparison writeup. Read this before touching solver internals. |
| `project_docs/` | Methodology notes (calibration) plus non-code reference artifacts (reports, spreadsheets, plots) from the physical bench/design side of the project. |
| `solver_design/` | Active and historical architecture design docs. **`FV_CORE_REWORK_PLAN.md` is the current plan of record** for the ongoing fluid-agnostic FV-core rework — read `solver_design/LLM_CONTEXT.md` for which other docs here are still current vs. superseded. |
| `EchTherm_correlation_docs/` | Extracted literature/reference-screenshot material backing the shell-and-tube (EchTherm-style) correlation set — Bell-Delaware shell-side and Ravigururajan-Bergles corrugated-tube correlations. |
| `reference/` | Large (~1 GB), gitignored tree of literature source PDFs/extracted-text backing the liquid/boiling coolant correlations. Do not enumerate file-by-file; see `reference/LLM_CONTEXT.md`. |
| `validation/` | Large, gitignored tree of generated validation-study run outputs (JSON/CSV/HTML/zip) written by scripts in `1Dmodel/validation/*.py`. Do not enumerate run folders; see `validation/LLM_CONTEXT.md`. Also holds `LIQUID_HEATED_CHANNEL_SOLVER_STATUS.md`, one of the files `CLAUDE.md` lists as canonical. |

## When to go where

- Need current solver architecture / pitfalls / wiring status → `context/`.
- Need calibration methodology → `project_docs/CALIBRATION_METHODOLOGY.md`.
- Need to know what's actively being built on the solver architecture, and
  what design docs are stale → `solver_design/` (start with
  `FV_CORE_REWORK_PLAN.md`'s status banner).
- Need a correlation's literature source for the shell-and-tube EchTherm
  comparison → `EchTherm_correlation_docs/`.
- Need a correlation's literature source for liquid/boiling coolant physics →
  `reference/` (do not walk it file-by-file; it's page-image/PDF extraction
  material, not code).
- Need evidence a validation run actually happened / its numeric result →
  `validation/` (do not walk generated run folders one by one; check the
  script that produced them in `1Dmodel/validation/` first).

## TODOs (only evidenced items)

None found stated at this level; see each subfolder's own file for its TODOs.

## Change history

No git history is meaningful here (repo has effectively one initial commit
plus a large amount of currently-uncommitted work). Dated notes live inside
individual docs, not here.
