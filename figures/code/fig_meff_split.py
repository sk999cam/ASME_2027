"""F15: regional split of ln(mu_eff / mu_eff,ref) into energy, length and alignment (paper Fig. F15, INTERIM).

Usage: python3 fig_meff_split.py
Output: /home/sk999/Documents/ASME_2027/figures/fig_meff_split.pdf

Three panels, one y scale each (the RANS terms are about 10 times the scale-resolving terms):
(a) RANS 50M against the 50M 3DNS, one bar per region;
(b) jet, near film (x/D 3-10), far film (x/D 10-25) and (c) diffuser, exit, approach layer: 20M LES (against the
50M 3DNS), 50M 3DNS and 100M 3DNS (against the 193M 3DNS), three bars per region.
Stacked bars: energy 1/2 ln k, length ln L_eff, alignment ln cos(theta). Diamond (right of the bar): total with its
full +-2 jackknife SE; the boxed number next to it is the total. The y limits of each panel contain every bar, marker
and error bar, so no error bar is clipped. Segment numbers are printed inside the segment and dropped when the segment
is shorter than the text.
Total marker convention: FILLED diamond = |total| >= 2 SE (resolved at about 95 %); HOLLOW diamond = |total| < 2 SE
(not resolved).
The script prints value and SE of every rung/region/term and lists every case with |value| < 2 SE.

INTERIM, the reference is not the 193M b010 for any rung. Data (analysis/ of /home/sk999/Documents/4dns_filmcooled/):
  RANS 50M vs 50M 3DNS : meff_ransp0_c50b/results_pp/regions_se.txt   (stat3 b006; point-wise regions)
  20M LES  vs 50M 3DNS : meff_c20bw_c50b/results_pp/regions_se.txt     (stat3 b006 of both; point-wise regions)
  50M 3DNS vs 193M     : meff_c50b_c193b/results/regions_se.txt        (b006 vs C193b b005, 5.0 FT; index boxes)
  100M 3DNS vs 193M    : meff_c100b_c193b/results/regions_se.txt       (same)
Work spec W1 (analysis/meff_<rung>_c193b10/results_pp/regions_se.txt, reference 193M b010) is used automatically per
rung if a file matching the glob in RUNGS exists (split_panels.pick); none existed on 2026-11-25. To point one rung at
a new file by hand, change its path in RUNGS.
The values cannot be chained between the rung groups (different point sets and selections). The density term is not
drawn (|density| <= 0.03), so the total can differ from the bar sum by that.
SE: delete-one-batch jackknife of both rungs in quadrature, (N-1)/N sum (t_i - mean)^2 (checked in
meff_c50b_c100b/errorbars.py: no unit, variance or factor error). It rests on 5-6 batches, so it is itself uncertain.
Style: set_mpl_params, Okabe-Ito term colours with a hatch per term (shared with fig_aeff_split.py), no title.
"""
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params
from split_panels import A, pick, read_se, report, panel, term_handles

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_meff_split.pdf'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
TERMS = [(r'energy $\frac{1}{2}\ln k$', '#56B4E9', ''), (r'length $\ln L_\mathrm{eff}$', '#E69F00', '///'),
         (r'alignment $\ln\cos\theta$', '#CC79A7', '...')]
COLS = [(c, h) for _, c, h in TERMS]
# (tick label, reference, file)
RUNGS = [('RANS 50M', 'vs 50M', pick(A + 'meff_ransp0_c50b/results_pp/regions_se.txt',
                                     'meff_*rans*50*_c193b10/results_pp/regions_se.txt')),
         ('20M LES', 'vs 50M', pick(A + 'meff_c20bw_c50b/results_pp/regions_se.txt',
                                    'meff_c20b*_c193b10/results_pp/regions_se.txt')),
         ('50M', 'vs 193M', pick(A + 'meff_c50b_c193b/results/regions_se.txt',
                                 'meff_c50b_c193b10/results_pp/regions_se.txt')),
         ('100M', 'vs 193M', pick(A + 'meff_c100b_c193b/results/regions_se.txt',
                                  'meff_c100b_c193b10/results_pp/regions_se.txt'))]
# (key in file, label)
REGIONS = [('jet, x/D 0–3', r'jet, $x/D$ 0–3'), ('film, x/D 3–10', r'near film, $x/D$ 3–10'),
           ('film, x/D 10–25', r'far film, $x/D$ 10–25'), ('hole, diffuser', 'diffuser'),
           ('hole exit', 'exit'), ('approach boundary layer', 'approach layer')]
YL = r'$\ln(\mu_\mathrm{eff}/\mu_\mathrm{eff,ref})$'


def entry(x, row):
    return (x, row[0], row[1], [row[2], row[4], row[6]])


data = [read_se(r[2]) for r in RUNGS]
report('F15 ln(mu_eff ratio)', [r[0] for r in RUNGS], data, [k for k, _ in REGIONS],
       [('total', 0), ('energy', 2), ('length', 4), ('align', 6)])

set_mpl_params(0.4, grid_on=True)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
fig, axs = plt.subplots(3, 1, figsize=(7.2, 8.4), gridspec_kw=dict(hspace=0.42))
fig.subplots_adjust(left=0.1, right=0.995, top=0.925, bottom=0.07)
axA, axB, axC = axs

# (a) RANS
panel(axA, [entry(g, data[0][k]) for g, (k, _) in enumerate(REGIONS)], COLS, YL, '(a) RANS 50M vs 50M 3DNS')
axA.set_xlim(-0.6, 6.0)
axA.set_xticks(range(6))
axA.set_xticklabels([n for _, n in REGIONS], fontsize=6.5)

# (b), (c): three rungs per region, bars at 4 g + j
for ax, gs, letter in ((axB, range(0, 3), '(b)'), (axC, range(3, 6), '(c)')):
    ent, ticks, labs = [], [], []
    for gi, g in enumerate(gs):
        k = REGIONS[g][0]
        for j in range(3):
            x = 4 * gi + j
            ent.append(entry(x, data[j + 1][k]))
            ticks.append(x)
            labs.append(RUNGS[j + 1][0] + '\n' + RUNGS[j + 1][1])
    panel(ax, ent, COLS, YL, letter)
    ax.set_xlim(-0.6, 11.0)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labs, fontsize=5.5)
    for gi, g in enumerate(gs):
        ax.annotate(REGIONS[g][1], (4 * gi + 1.15, 0), xycoords=('data', 'axes fraction'), textcoords='offset points',
                    xytext=(0, -26), ha='center', va='top', fontsize=7)
        if gi:
            ax.axvline(4 * gi - 0.45, color='0.6', lw=0.6, ls='--')

handles = term_handles(TERMS)
handles.append(Line2D([], [], marker='D', color='k', ls='', label=r'total, $\pm2$ SE'))
handles.append(Line2D([], [], marker='D', color='k', mfc='white', ls='', label=r'total, $|\mathrm{total}|<2$ SE'))
fig.legend(handles=handles, loc='lower center', ncol=3, bbox_to_anchor=(0.55, 0.935))
fig.savefig(OUT, bbox_inches='tight')
print('wrote', OUT)
