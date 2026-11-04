# Log

Newest entry first. One entry per working session: date, what changed, files touched.
For each added figure, record the source script and data in 4dns_filmcooled.

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
