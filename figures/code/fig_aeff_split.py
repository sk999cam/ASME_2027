"""F16: regional split of ln(alpha_eff / alpha_eff,ref) and ln(Pr_t,eff / Pr_t,eff,ref) (paper Fig. F16, INTERIM).

Usage: python3 fig_aeff_split.py
Output: /home/sk999/Documents/ASME_2027/figures/fig_aeff_split.pdf

Four panels, one y scale each (RANS terms are 10-30 times the LES terms), one bar per region:
(a) RANS 50M and (b) 20M LES: ln(alpha_eff/alpha_eff,ref), stacked energy 1/2 ln k, thermal length ln L_T and
    alignment ln cos(phi) terms. The thermal mixing length l_T and the correlation rho_uT are NOT stored separately:
    the analysis files give only the combined L_T = alpha_eff / (rho k^1/2 cos phi), so L_T is drawn and labelled as such.
(c) RANS 50M and (d) 20M LES: ln(Pr_t,eff/Pr_t,eff,ref) = ln(L_eff/L_T) [length ratio] + ln(cos theta/cos phi)
    [alignment ratio]. Under each region name of (c), the reference Pr_t,eff = 0.9 exp(-total_RANS): the RANS has
    Pr_t = 0.9 exactly (cos theta = cos phi = 1), so the reference value follows from the RANS ratio. It is the
    volume-averaged-ln (geometric) value, not the mean of Pr_t,eff.
Diamond (right of the bar): total with its full +-2 jackknife SE; the boxed number next to it is the total. The y
limits of each panel contain every bar, marker and error bar, so no error bar is clipped. Segment numbers are dropped
when the segment is shorter than the text (all numbers are printed on stdout).
Total marker convention: FILLED diamond = |total| >= 2 SE (resolved at about 95 %); HOLLOW diamond = |total| < 2 SE.
The script prints value and SE of every rung/region/term and lists every case with |value| < 2 SE.

INTERIM: both rungs are against the 50M 3DNS, not the 193M. No alpha_eff / Pr_t,eff comparison with the 193M exists
(only mu_eff), so the 50M and 100M rungs of F15 have no bars here. Chaining through the 50M is not valid. Data, in
/home/sk999/Documents/4dns_filmcooled/analysis/ (point-wise regions, stat3 b006 of both rungs):
  meff_ransp0_c50b/results_pp/regions_alpha_se.txt, regions_prt_se.txt   (RANS p0dns 50M vs C50b)
  meff_c20bw_c50b/results_pp/regions_alpha_se.txt, regions_prt_se.txt    (C20bW vs C50b)
Work spec W1 (analysis/meff_<rung>_c193b10/results_pp/regions_{alpha,prt}_se.txt, reference 193M b010) is used
automatically if a file matching the glob in RUNGS exists (split_panels.pick); none existed on 2026-11-25. To point a
rung at a new folder by hand, change its entry in RUNGS (folder holding both files).
The approach layer is threshold-dependent (volume used 0.11-0.13 of the region; FINDINGS C-threshold, write-up
03_differences) and is marked with *. The density term (|d| <= 0.03) is not drawn; the total is the file total.
SE: delete-one-batch jackknife of both rungs in quadrature (see fig_meff_split.py; no unit or variance error).
Style: set_mpl_params, same term colours and hatches as fig_meff_split.py, no title.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params
from split_panels import A, pick, read_se, report, panel, term_handles

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_aeff_split.pdf'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
C_EN, C_LEN, C_AL = '#56B4E9', '#E69F00', '#CC79A7'
H_EN, H_LEN, H_AL = '', '///', '...'
PR_RANS = 0.9


def folder(interim, w1_glob):
    """W1 folder if its alpha table exists, else the interim folder."""
    f = pick(interim + 'regions_alpha_se.txt', w1_glob + '/regions_alpha_se.txt')
    return os.path.dirname(f) + '/'


RUNGS = [('RANS 50M', folder(A + 'meff_ransp0_c50b/results_pp/', 'meff_*rans*50*_c193b10/results_pp')),
         ('20M LES', folder(A + 'meff_c20bw_c50b/results_pp/', 'meff_c20b*_c193b10/results_pp'))]
REGIONS = [('jet, x/D 0–3', 'jet'), ('film, x/D 3–10', 'near\nfilm'), ('film, x/D 10–25', 'far\nfilm'),
           ('hole, diffuser', 'diffuser'), ('hole exit', 'exit'), ('approach boundary layer', 'approach\nlayer*')]
KEYS = [k for k, _ in REGIONS]

alpha = [read_se(r[1] + 'regions_alpha_se.txt') for r in RUNGS]
prt = [read_se(r[1] + 'regions_prt_se.txt') for r in RUNGS]
names = [r[0] for r in RUNGS]
report('F16 ln(alpha_eff ratio)', names, alpha, KEYS, [('total', 0), ('energy', 2), ('L_T', 4), ('align', 6)])
report('F16 ln(Pr_t,eff ratio)', names, prt, KEYS, [('total', 0), ('length', 2), ('align', 4)])

set_mpl_params(0.4, grid_on=True)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.0), gridspec_kw=dict(hspace=0.5, wspace=0.3))
fig.subplots_adjust(left=0.1, right=0.995, top=0.88, bottom=0.1)
YA = r'$\ln(\alpha_\mathrm{eff}/\alpha_\mathrm{eff,ref})$'
YP = r'$\ln(\mathrm{Pr}_{t,\mathrm{eff}}/\mathrm{Pr}_{t,\mathrm{eff,ref}})$'
ref = [PR_RANS * np.exp(-prt[0][k][0]) for k in KEYS]

for j in range(2):
    ax = axs[0, j]
    panel(ax, [(g, alpha[j][k][0], alpha[j][k][1], [alpha[j][k][2], alpha[j][k][4], alpha[j][k][6]])
               for g, k in enumerate(KEYS)],
          [(C_EN, H_EN), (C_LEN, H_LEN), (C_AL, H_AL)], YA, f'({"ab"[j]}) {names[j]}')
    ax = axs[1, j]
    panel(ax, [(g, prt[j][k][0], prt[j][k][1], [prt[j][k][2], prt[j][k][4]]) for g, k in enumerate(KEYS)],
          [(C_LEN, H_LEN), (C_AL, H_AL)], YP, f'({"cd"[j]}) {names[j]}')
for ax in axs.ravel():
    ax.set_xlim(-0.6, 6.4)
    ax.set_xticks(range(6))
for j in range(2):
    axs[0, j].set_xticklabels([n for _, n in REGIONS], fontsize=5.5)
    axs[1, j].set_xticklabels([n + ('\n' + rf'ref {r:.2f}' if j == 0 else '') for (_, n), r in zip(REGIONS, ref)],
                              fontsize=5.5)

handles = term_handles([(r'energy $\frac{1}{2}\ln k$', C_EN, H_EN),
                        (r'length: (a, b) $\ln L_T$, (c, d) $\ln(L_\mathrm{eff}/L_T)$', C_LEN, H_LEN),
                        (r'alignment: (a, b) $\ln\cos\phi$, (c, d) $\ln(\cos\theta/\cos\phi)$', C_AL, H_AL)])
handles.append(Line2D([], [], marker='D', color='k', ls='', label=r'total, $\pm2$ SE'))
handles.append(Line2D([], [], marker='D', color='k', mfc='white', ls='', label=r'total, $|\mathrm{total}|<2$ SE'))
handles.append(Line2D([], [], ls='', label='* threshold-dependent region; ref: reference $\\mathrm{Pr}_{t,\\mathrm{eff}}$'))
fig.legend(handles=handles, loc='lower center', ncol=2, bbox_to_anchor=(0.55, 0.89))
fig.savefig(OUT, bbox_inches='tight')
print('wrote', OUT)
