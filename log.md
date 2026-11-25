# Log

Newest entry first. One entry per working session: date, what changed, files touched.
For each added figure, record the source script and data in 4dns_filmcooled.

## 2026-11-25 (3): figure fixes, Pr_t story set to a film/jet split

- Redrawn at 193M b010: F5, F6 (now hole sections at s/D −2.3, −1.1, −0.5), F13, F14, F17, F18, F19.
  Threshold masks removed (F14, F17, F18, F19); F17 L_eff on a log scale; F14 has a dashed threshold line.
- F15, F16: labels placed on their bars, full ±2 SE bars, hollow marker when |total| < 2 SE
  (`figures/code/split_panels.py`). No SE bug.
- 193M b010 is the cumulative 10 FT mean (mean_time 0.0359 s = 400 D/U). Near-wall film layer (z/D < 0.3)
  not converged in α_eff: sign changes at 31 % of points between b006 and b010.
- Owner decision: heat-flux error = momentum error + Pr_t error (E12); classic Pr_t (E10b) shown in the film
  against Kohli and Bogard; no Pr_t quoted for the jet, where a scalar Pr_t does not fit; regional values from
  a gradient-weighted least-squares fit (E15). Plan 4.1 (iv), 4.2 (i), (ii), F14, F16, F20 rewritten.
- e9: W1 amended (E12 terms, least-squares estimator, R²); asked what the split means for the corrected RANS
  runs. Pr_t class analysis running (`fig_prt_classes.py`, `fig_prt_classic.py`).

## 2026-11-25 (2): figure feedback round 1

- Comments taken from `planning/feedback/writeup_plan_annotated_2026-11-25.pdf` (7 notes).
- Reference changed to the 193M at b010 for every figure and number (user). Lower rungs stay at b006. No empty
  panels and no stand-in rungs (user). Brief updated: `planning/figure_brief_common.md`.
- Plan changes:
  - F7 gets a full explanation of panels (a)–(c), including the 95 % coherence level.
  - The RANS "low-γ sublayer" remark is removed (a viscous sublayer exists in the DNS too).
  - New F9 (`fig_spots.py`): space-time diagram and snapshot sequence as the evidence for spot growth and merger.
    Old F9–F19 renumbered to F10–F20.
  - New bullet on the mid-pitch η gap against the experiment. Hypotheses: the tripped turbulent approach layer of
    the experiment (H 1.45, Re_θ 670 at x/D −4.7, Schroeder & Thole 2014) against the laminar 3DNS layer; and the
    distance transition needs. Test: γ_mid overlaid on η_mid (F12 (c)).
- Calculations dispatched:
  - Three Sonnet agents:
    - hole and η-wall at b010: F5, F6 (centre plane replaced by a close-up), F13;
    - plane maps and budget lines at b010: F14, F17, F18, F19 (unique panel letters);
    - γ for the 193M and the 20M LES, F8–F12 at b010, and the new F9.
  - Work spec to e9 (`planning/work_spec_e9_2026-11-25.md`):
    - W1: regional splits against the 193M b010 for F15 and F16;
    - W2: directional diffusivities for the 193M and 100M, for F20;
    - W3: bubble motion on the 193M snapshots, for F7.
    Due 3 Dec.

## 2026-11-25: plan feedback from the PDF; Results 4.1 and 4.2 planned

- User comments read from the annotated `planning/writeup_plan.pdf` (pypdf); the annotated copy is kept as
  `planning/feedback/writeup_plan_annotated_2026-11-24.pdf`, since a rebuild overwrites the PDF.
  - P1: reference for RANS as the design tool: Milani et al. 2018, 2020.
  - P4: gap rewritten. The reference is taken as resolved (3.3); the missing piece is the split of the RANS error
    by region and cause and its link to the eta-bar error.
  - P5: rewritten from the literature. Fawcett (2011): no coherent unsteadiness in the shaped-hole jet at
    BR 0.8-1.0 (DR 1, VR near 1). The present VR is about 1. The gap is the share of the mixing physics carried
    by each mechanism, not the resolved unsteadiness.
  - P6: T1 (prior studies) dropped. Replaced by one paragraph grouping the studies by resolution and modelling.
    The rungs table is now T1 (`\setcounter{table}{0}`). P8 justification chain updated.
- Results planned (bullets only): 4.1 Coolant Flow (hole chain, mainstream chain, effectiveness, the
  mu_eff/alpha_eff/Pr_t,eff frame) and 4.2 Breakdown of Model Differences (three-factor split, log budget,
  region by region, literature, remedies). 4.3 A Posteriori Tests is a placeholder. New IDs: F5-F19, E9-E16.
  Numbers from `4dns_filmcooled/analysis/FINDINGS.md` (H, S, LS, SU, V, C, R series) and
  `write_up/results/02_coolant_flow.tex`, `03_differences.tex`. Rungs: the six of the Methods; no corrected RANS.
  Reference 193M b006. The plan PDF is 14 pages.

## 2026-11-24: write-up plan, Part 1 (Methods), revision 2

- `planning/writeup_plan.tex/.pdf`: section-by-section plan of the Methods. Three subsections: 2.1
  Configuration, Computational Methods and Boundary Conditions (with the analysis regions), 2.2 Mesh Family,
  2.3 Resolution. Resolution stated as Delta_mesh/eta_193 with eta from the 193M 3DNS batch b010 (10 FT).
- `planning/plan_figures.py` writes `planning/figures/fig_structure.pdf`, `fig_config_regions.pdf` (draft of
  paper Fig. 1: domain, BCs and regions), `fig_ladder_map.pdf`. Region bounds from
  `4dns_filmcooled/dns/python/resolution_regions.py`; hole outline from `write_up/methods/01_configuration.tex`.
- Open decisions are listed in section 4 of the plan.
- Revision 3 (user): ladder map and cell-size figure dropped; T1 is now a rungs-and-meshes table with s;
  5M is a family member (s ~0.30 from 5.33M cells, counts to be read from `rans_5M_collar22.cgns`); the
  resolution figure and table cover 5M-193M; downstream regions renamed near film / far film.
- Figures remade to the CLAUDE.md rules and saved in `figures/` (plan reads them via `\graphicspath{{../figures/}}`;
  `planning/figures/` removed):
  - `fig_config_regions.pdf` from `4dns_filmcooled/write_up/methods/python/fig_config_regions.py` (new script).
  - `fig_mesh_classes.pdf` from `plot_mesh_classes.py`, mesh `/mnt/Elements/DNS/meshes/rans_20M_Lm25b.cgns`.
  - `fig_resolution.pdf` from `fig_resolution.py`: interim, each run's own eta; batches C20bW b006, C50b b006,
    C100b b006, C193b b010. Panels (user): diffuser, exit, near film (x/D 3-10), far film (x/D 10-25),
    whole domain; metering left out.
- `figures/` holds only candidate paper figures; each has its script in `figures/code/` with the same name
  (`fig_config_regions.py`, `fig_mesh_classes.py` = copy of `plot_mesh_classes.py`, `fig_resolution.py`).
  Scripts use absolute paths and run from there without arguments; all three were rerun from `figures/code/`.
- Plan: Methods equations E1-E6 added (Re_D, BR/DR/VR/I, effectiveness, N(s), xi_m, Delta_mesh/eta_193) with
  purpose of each; open decision on the eta symbol clash (effectiveness vs Kolmogorov length).
  `planning/build_plan.sh` rebuilds the PDF (`--figs` also redraws figures).
- Plan revision 4: all decisions taken (listed in section 4 of the plan); no pending items. Headings 2.1
  Simulation Configuration and Computational Methods, 2.2 Mesh Family and Coarsening Method, 2.3 Achieved Mesh
  Resolution; eta_K for the Kolmogorov length. T1 cells from the CGNS files (5M 5.77M ... 193M 192.3M); BR from
  `analysis/hole_discharge/README.md` (common p0); 5M s = 0.30 (point-weighted block average). F3 axis labels
  (user edit) kept, eta -> eta_K; original script synced to figures/code copy.
- Boundary-condition table dropped (user); conditions on F1(a) and in one paragraph. Resolution table renumbered T3.
- Boundary conditions stated once in text, the same for every rung (user). BR matching removed from the Methods:
  introduced in the Results with the variety in diffuser-bubble behaviour (user).
- `writeup_plan` made the plan for the whole paper (user): figure/table/equation lists at the start; decisions
  section and Methods structure section removed (decisions kept in chat and here); `plan_figures.py` and
  `fig_plan_methods_structure.pdf` deleted; `build_plan.sh` updated.
  Decisions so far: machine/ranks not reported; F1 full width; wall units one sentence; Delta mapped onto 193M;
  eta effectiveness / eta_K Kolmogorov; headings 2.1 Simulation Configuration and Computational Methods, 2.2 Mesh
  Family and Coarsening Method, 2.3 Achieved Mesh Resolution; 5M s 0.30; BR common p0, BR matching in Results;
  no BC table, BCs in text.
- T1 filled (status draft): 12 prior studies plus the present study, read from the PDFs in
  `4dns_filmcooled/write_up/literature/papers/` (pdftotext). Columns: hole, method, cells, approach layer, BR,
  resolution check, closure assessed. Mesh sizes marked "p" are grid points. Not found in the paper: Yu and
  Yavuzkurt (2018) approach-layer state; Nunno et al. (2022) and Tyagi and Acharya (2003) resolution check.
- Introduction decisions (user): F1 reused from Kar et al. 2026 Fig. 1 (vector crop by
  `figures/code/fig_transition_map.py`, not redrawn); T1 kept; IDDES not a rung; aims cover performance,
  resolution and closure; RANS failure leads; E1/E2 in the Introduction; ~2 pages; aims in a separate section
  2 Research Aims (as in the GPPS paper), so Methods is Section 3 (3.1-3.3). New P5 on unsteady flow structures
  and the hole flow; paragraphs now P1-P8.
- Introduction planned (bullets only): P1 film cooling and shaped holes, P2 7-7-7 baseline, P3 flat plate at ZPG,
  P4 RANS errors (E1 Boussinesq, E2 GDH), P5 scale-resolving simulations (T1 table of prior studies), P6 laminar
  approach layer (F1 transition map, reproduced from Clark et al. 2024; source data not on disk), P7 aims and
  justification chain into Methods. IDs renumbered in paper order: Methods figures F2-F4, table T2, equations E3-E8.
- Class table (T2) and resolution table (T3) dropped (user); their content goes in the text. Methods keeps one
  table, T1 (rungs and meshes).
- Plan-only diagram `fig_plan_methods_structure.pdf` moved to `planning/` (later deleted) (drawn by `planning/plan_figures.py`).

## 2026-11-04: kar2026 PDF link

- `references.bib`, `kar2026`: added `url` = https://gpps.global/wp-content/uploads/2026/09/GPPS-TC-2026_005.pdf
  (HTTP 200, application/pdf, 2 973 676 bytes, same size as the local copy). The reference list prints
  both the DOI and the URL; the URL gives a working link while the DOI is unregistered.

## 2026-11-04: reference group added to the bibliography

- `4dns_filmcooled/write_up/literature/refs.json`: new group `reference` with the 2 papers in
  `papers/reference/`: 10.33737/gpps26-tc-005 (Kar et al. 2026, GPPS Montreal26) and
  10.1115/GT2024-124918 (Clark et al. 2024). That repo's 04:05 sync commits the change.
- `manuscript/references.bib`: `kar2026` and `clark2024` written by hand in asmeconf `@inproceedings`
  format (paper number, venue, event date). doi.org returned "DOI Not Found" for 10.33737/gpps26-tc-005
  on 2026-11-04, so its link will not resolve until GPPS registers it. The PDF gives the paper number
  as GPPS-TC-2026-0005.
- `main.tex` abstract: "(Kar et al., 2026)" replaced by `\cite{kar2026}`.
- BibTeX with `asmeconf.bst`: 88 items, 0 errors.

## 2026-11-04: abstract and paper number

- `main.tex`: paper number set to GT2027-202699.
- `main.tex`: abstract inserted from `write_up/framework/abstract.txt` (2026-10-15). Wording unchanged;
  LaTeX markup only (`\Reth`, en dashes, `\%`, math for $k$--$\omega$ and $x/D$). Schroeder and Thole
  (2014) replaced by `\cite{schroeder2014}`. Kar et al. (2026) kept as plain text: no bib entry yet.
  Length 401 words against the ASME limit of 200.

## 2026-11-04: sync order and instructions

- `git_sync.sh`: order changed to commit → pull → push (was pull → commit → push). An uncommitted local
  edit to a file that also changed on Overleaf made the old pull fail every night.
- `plan.md` §3: sync instructions for Overleaf ↔ GitHub ↔ local.
- First Overleaf push (a0dafbd, affiliation edit in `main.tex`) also reset `git_sync.sh` from mode 755
  to 644. Cron entry changed to `/bin/bash /home/sk999/Documents/ASME_2027/git_sync.sh` so the mode no
  longer matters. Mode restored to 755 in 404d5e0.

## 2026-11-04: repository set up

- Created `plan.md`, `log.md`, `.gitignore`.
- `manuscript/main.tex`: skeleton built on asmeconf v1.47 with GT2027 header placeholders, working title,
  keywords from `write_up/framework/abstract.txt`, and section outline from
  `write_up/framework/paper_structure.tex`. Class options add `lineno` for drafting. Language options
  and sample content from the template removed. `\graphicspath` points to `../figures/`.
- `manuscript/asmeconf.cls`, `manuscript/asmeconf.bst`: copied unchanged from `template/`.
- `manuscript/references.bib`: 86 entries. 83 from doi.org (Crossref BibTeX) for the DOIs in
  `4dns_filmcooled/write_up/literature/refs.json` (background 40, methods 37, conference versions 6).
  3 arXiv preprints (2105.06225, 2502.19609, 2503.03844) entered by hand in the asmeconf `eprint` format.
  Unicode dashes converted to `--`/`---`, `url`, `ISSN` and `month` fields dropped.
  BibTeX with `asmeconf.bst`: 86 items, 0 errors, 38 warnings, all missing `pages` (ASME journal
  articles have no page range).
- `template/` git-ignored and left unchanged.
- Local compile of `main.tex` fails: TeX Live 2018 here lacks the expl3 commands asmeconf v1.47 uses.
  Overleaf compiles it. BibTeX runs locally.
- `git_sync.sh`: nightly pull, commit, push at 04:10 (cron), log in `git_sync.log`. Origin switched from
  HTTPS to SSH so cron runs without a credential prompt.
