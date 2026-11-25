# Common brief for the figure agents (GT2027 paper, Results 4.1 and 4.2)

You are drawing candidate figures for an ASME GT2027 paper about a 7-7-7 laidback fan-shaped film-cooling hole.
The configuration is DR 1.5 and BR about 1.5, with a laminar approach layer under freestream turbulence, on a flat
plate at zero pressure gradient.

The plan is `/home/sk999/Documents/ASME_2027/planning/writeup_plan.tex`, Section 4. Read the entries for your
figures before you start. The figure specs below are the minimum.

## Rungs (the only cases to show)

| label in figures | case | data |
|---|---|---|
| RANS 5M | `starccm/rans_5M_Lm25_BR15_FST_ZPG_p0dns` (SST + γ–Re_θ, steady) | STAR exports / post files in that case folder |
| RANS 50M | `starccm/rans_50M_Lm25b_BR15_FST_ZPG_p0dns` | as above (cell exports exist, see `analysis/plume`, `analysis/layer_state`) |
| 20M LES | C20bW | `/mnt/Tank/DNS/film_cooling/C20bW_BR15/pictures_job3/stat3_b006` |
| 50M 3DNS | C50b | `/mnt/Tank/DNS/film_cooling/C50b_BR15/pictures_job9/stat3_b006` |
| 100M 3DNS | C100b | `/mnt/Elements/DNS/film_cooling/C100b_BR15/pictures_job9/stat3_b006` |
| 193M 3DNS (reference) | C193b | **batch b010 (10 FT): `/mnt/Elements/DNS/film_cooling/C193b_BR15/pictures_job18/stat3_b010`** |

- Paths starting with `starccm/`, `analysis/` or `dns/` are relative to `/home/sk999/Documents/4dns_filmcooled/`.
- Do not show the IDDES, the 4040 Pa RANS, the BR-matched RANS, or any QCR / Pr_t / substitution RANS.
- The reference is ALWAYS the 193M 3DNS at batch b010 (user, 2026-11-25). Never use 193M b006 or b005. The lower
  rungs keep their last batch, b006. Time series (samplers, snapshots, monitor) of the 193M use the b010 window:
  take its start and end time from `stat3_batches.txt` in the 193M run folder.
- No empty panels and no stand-in rungs (user, 2026-11-25). If a quantity has not been computed for a rung,
  compute it: that is part of your task. If it cannot be computed (no raw data exists), stop and report this.
  Do not draw another rung in its place. Never invent, interpolate across rungs or fake data.
- Regions are named metering, diffuser, exit, approach layer, jet (x/D 0–3), near film (x/D 3–10), far film
  (x/D 10–25) and freestream. Write "near film" and "far film", never "wake". Analysis files often use
  `wake_3_10` / `wake_10_25` for these.

## Where to find things

- Findings with numbers and sources: `analysis/FINDINGS.md`. Per-analysis READMEs: `analysis/<folder>/README.md`.
  Existing plotting scripts are in the same folders.
- Results drafts that already use these data: `write_up/results/02_coolant_flow.tex` and `03_differences.tex`.
- 3DNS readers and tools: `dns/python/` (e.g. `plot_stat3.py`, `plot_stat3_planes_pf.py`, `tke_budget_pf.py`,
  `length_scales_pf.py`, `resolution_regions.py`), `analysis/meff_c50b_c100b/dns_fields.py`, `compare.py`.
- Reuse existing code by importing it or copying it into your script. Do not edit any file in `4dns_filmcooled`.

## Figure rules (from the project CLAUDE.md; these are mandatory)

1. Before any figure is created, call `set_mpl_params(sf, grid_on=...)` from
   `/home/sk999/Documents/3DNS_FP/functions.py` (`sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')`).
   Do not override its rcParams (fonts, font sizes, line widths, marker sizes, grid, legend, spines). Change
   `sf` or `figsize` instead.
2. Colours: Okabe-Ito only, set explicitly after `set_mpl_params`:
   `plt.rcParams['axes.prop_cycle'] = cycler(color=['#000000','#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7'])`.
   Distinguish series by line style as well (`-`, `--`, `-.`, `:`). Markers only `o s ^ v D x +`.
   - Use one fixed colour and style per rung in every figure:

     | rung | colour | style |
     |---|---|---|
     | 193M | `#000000` | `-` |
     | 100M | `#0072B2` | `--` |
     | 50M | `#D55E00` | `:` |
     | 20M LES | `#E69F00` | `-.` |
     | RANS 50M | `#009E73` | `-` |
     | RANS 5M | `#CC79A7` | `--` |
     | experiment | black markers `o` |
3. No titles: no `plt.title`, `ax.set_title` or `fig.suptitle`. Label panels (a), (b), ... inside the axes.
4. Colour maps:
   - sequential data: `inferno`;
   - `jet` only for temperature and temperature-related quantities (η, Pr_t is not temperature);
   - signed data and log ratios: `bwr` with limits symmetric about zero;
   - no other colour maps.
5. Formats: `.pdf` for line plots and diagrams, `.png` for contour and field plots (dpi ≥ 300).
6. Output: save directly to `/home/sk999/Documents/ASME_2027/figures/<name>.<ext>`. The script is
   `/home/sk999/Documents/ASME_2027/figures/code/<name>.py`, with the same name, one script per figure.
   - The script runs with `python3 <name>.py` and no arguments, using absolute paths.
   - Its docstring states: what the figure shows, the data files and batches, and anything interim.
7. Axis labels with units or normalisation, e.g. `$x/D$`, `$k/U_\infty^2$`, `$s/D$`. Legends without boxes are
   fine.

## Practicalities

- Shared machine: 125 GB RAM and 32 cores, with several agents running at once.
  - Use at most about 20 GB RAM and 4 processes.
  - Read 193M / 100M data block by block.
  - `/mnt/Elements` reads at about 10–34 MB/s: read only the stat3 columns and blocks you need.
- Cache intermediate arrays (small `.npz`, at most a few GB) in
  `/mnt/Elements/DNS/analysis/asme2027_figs/<name>/`. Make the script reuse the cache if it exists. Never
  write to `/mnt/Tank` (98 % full).
- Write files with your file-writing tools. Use Python only to run scripts.
- After saving each figure, open the image (Read tool) and check it: labels readable, nothing clipped, no title,
  colours as above. Fix it and redraw if needed.

## Report back (final message, terse)

For each figure:
- path;
- data used (cases, batches, files);
- status (complete / interim, and why);
- 2–4 key numbers the figure shows;
- anything in the spec you could not do, and why.
