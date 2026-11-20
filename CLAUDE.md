# Rules for the ASME 2027 paper

## Writing

- Any text written for the paper or about it (manuscript sections, captions, responses to reviewers) follows the
  author's style. Read `/home/sk999/Documents/sample_writing/WRITING_STYLE.md` before writing.
- Passive voice always. Never use "we", "I", "our", "us" or "my".
- Exponents are written as `$3\times10^{-6}$`. No e-notation.
- Do not re-read the PDFs in `/home/sk999/Documents/sample_writing/` unless the style guide does not cover a case.
- Keep the author's voice. Fix the faults listed in section 3 of the style guide while writing; do not rewrite the
  voice towards a different style.

# Figure rules

These rules apply to every figure produced for the paper.

## 0. Figure pipeline

- The pipeline is fixed: `/home/sk999/Documents/4dns_filmcooled` (data and plotting scripts) ->
  `/home/sk999/Documents/ASME_2027/figures` (saved figures) -> the `.tex` manuscript.
- Every figure used in the paper is saved directly into `/home/sk999/Documents/ASME_2027/figures`. Plotting scripts
  set their output path to this folder.
- Do not save paper figures anywhere in `4dns_filmcooled`, and do not leave copies scattered across it.
- The `.tex` file includes figures only from `figures/`.
- Format: `.pdf` for line plots and diagrams, `.png` for contour and field plots, `.jpeg` only for photographs.

## 1. Plot style: always call `set_mpl_params`

- Every plotting script calls `set_mpl_params(sf, grid_on=...)` from
  `/home/sk999/Documents/3DNS_FP/functions.py` before creating any figure.
- Do not override its rcParams (fonts, font sizes, line widths, marker sizes, grid, legend, spines) in individual
  scripts. If a figure needs a different size, change `sf` or `figsize`, not the individual params.

## 2. Colour-blind-friendly line colours

- Line and marker colours come from a colour-blind-safe palette. Use the Okabe-Ito set:
  `#000000, #E69F00, #56B4E9, #009E73, #F0E442, #0072B2, #D55E00, #CC79A7`.
- `set_mpl_params` loads the `seaborn-v0_8-bright` style, whose default colour cycle is not colour-blind safe. Set
  the colour cycle explicitly after calling it, e.g.
  `plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)`.
- Distinguish series with line styles as well as colour: solid, dashed (`--`), dash-dot (`-.`), dotted (`:`).
- Markers are allowed. Keep to standard shapes: `o`, `s`, `^`, `v`, `D`, `x`, `+`. No unusual or custom markers.

## 3. No titles

- No figure going into the paper has a title. Do not call `plt.title`, `ax.set_title` or `fig.suptitle`. The caption
  in the manuscript describes the figure.

## 4. Colour maps

- One-sided (sequential) data: `inferno` or `jet` only.
- Diverging data (signed about zero): `bwr` only, with limits symmetric about zero.
- `jet` is reserved for temperature and temperature-related quantities. All other sequential data use `inferno`.
- No other colour maps.
