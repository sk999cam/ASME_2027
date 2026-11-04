# Plan: ASME Turbo Expo 2027 paper

Working title: *Resolution Ladder for the 7-7-7 Laidback Fan-Shaped Film-Cooling Hole*.
The paper reports the 4dns_filmcooled campaign: RANS, URANS, LES and DNS of the 7-7-7 hole at
DR 1.5, BR 1.5, laminar approach boundary layer.

Record decisions and status changes here. Record dated work in `log.md`.

## 1. Workflow

1. Draft the conceptual version of each section in `../4dns_filmcooled/write_up/`
   (`framework/`, `methods/`, `results/`, `literature/`).
2. Soham writes the manuscript text in Overleaf from those drafts. The manuscript is Soham's own writing.
3. Claude works on request only:
   - rewrites a named section for clarity;
   - generates figures from 4dns_filmcooled data and inserts them;
   - reviews Soham's text and gives an overview;
   - fixes LaTeX formatting (asmeconf class, equations, tables, references).
4. Overleaf syncs with GitHub (`sk999cam/ASME_2027`). `git_sync.sh` syncs this clone with GitHub nightly.

## 2. Repository layout

| path | content | in Overleaf |
|---|---|---|
| `manuscript/main.tex` | the paper (asmeconf v1.47) | yes, main document |
| `manuscript/references.bib` | bibliography, `asmeconf.bst` style | yes |
| `manuscript/asmeconf.cls`, `asmeconf.bst` | class and style, copied unchanged from `template/` | yes |
| `figures/` | final figures (PDF for line plots, PNG at 300 dpi or more for contours) | yes |
| `plan.md`, `log.md` | this plan and the work log | yes (harmless) |
| `git_sync.sh` | nightly sync script | yes (harmless) |
| `template/` | asmeconf distribution, reference only, never edited | no (git-ignored) |

Rules:
- Copy anything needed from `template/` into `manuscript/` before changing it. `template/` stays as downloaded.
- Raw data and post-processing stay in `4dns_filmcooled`. Only finished figures come here.
- Figure scripts stay in `4dns_filmcooled` next to the data they read. Each figure's source script
  goes in `log.md` when the figure is added.
- Keep files under 50 MB (Overleaf upload limit). The sync script skips larger files.

## 3. Overleaf connection

1. Overleaf: New Project → Import from GitHub → `sk999cam/ASME_2027`.
2. Menu → Main document: `manuscript/main.tex`. Compiler: pdfLaTeX, TeX Live 2025 or later
   (asmeconf v1.47 does not compile on the local TeX Live 2018).
3. Edits made on Overleaf reach GitHub only through Overleaf's Menu → GitHub → Push. Edits pushed to
   GitHub reach Overleaf only through Menu → GitHub → Pull. Push from Overleaf before the 04:10 sync
   and pull after it to avoid merge conflicts.

## 4. Paper structure and status

Outline from `4dns_filmcooled/write_up/framework/paper_structure.tex`.

| section | conceptual draft (write_up) | manuscript |
|---|---|---|
| Abstract | `framework/abstract.txt` | not started |
| 1 Introduction | `literature/`, `framework/mechanisms_from_literature.md` | not started |
| 2.1 Configuration and operating point | `methods/01_configuration.tex` | not started |
| 2.2 Simulation ladder | `methods/02_simulation_ladder.tex` | not started |
| 2.3 Validation | — | not started |
| 2.4 Statistics and sampling | `methods/03_mesh_resolution.tex` | not started |
| 2.5 Analysis method | `methods/04_analysis.tex`, `framework/paper_maths.tex` | not started |
| 3.1 Performance across the ladder | `results/` | not started |
| 3.2 Departures of modelled from resolved turbulence | `results/`, `framework/closure-error-analysis-plan.md` | not started |
| 3.3 Energy, length scale and alignment | `framework/eddy_viscosity_link.md`, `framework/mu_eff_origins.tex` | not started |
| 3.4 Flow structures and mechanisms | `results/` | not started |
| 3.5 Departures that matter for effectiveness | `results/` | not started |
| 3.6 Resolution requirements and cost | `results/` | not started |
| 4 Conclusions | — | not started |

## 5. Deadlines

TODO: enter the Turbo Expo 2027 dates (abstract, draft paper, review return, final paper) from the call.

| item | date |
|---|---|
| abstract | |
| full draft for review | |
| final paper | |

## 6. Open items

- Conference header fields in `main.tex` (date, city, paper number) are placeholders.
- Co-authors, affiliations and corresponding e-mail in `main.tex` are placeholders.
- `references.bib` titles come from Crossref as published. ASME wants title case: check each cited entry.
- Remove the `lineno` and `colorlinks` class options before final submission.
