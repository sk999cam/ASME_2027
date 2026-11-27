"""F15: regional split of ln(mu_R,RANS / mu_R,193), gradient-weighted least squares (mode C), two panels.

Usage: python3 parse_meffR_log.py /home/sk999/Documents/ASME_2027/planning/maths_audit/estimator/estimator_modeC_x20.log
       (once, writes data/meffR_modeC.txt and data/meffR_se.txt)
       python3 fig_meff_split.py
Output: /home/sk999/Documents/ASME_2027/figures/fig_meff_split.pdf

Data: data/meffR_modeC.txt (from planning/maths_audit/estimator/estimator_modeC.log) and data/meffR_se.txt.
The SE of the total is read from data/meffR_se.txt (columns region, SE; the line '# kind: ...' names the SE kind and
goes into the legend), so a time-batch SE can replace the spatial-block one by overwriting that file.
mu_R = gradient-weighted least squares, wall nodes excluded, z/D >= 0.3. Regions: jet x/D 0-3, near film 3-10, far
film 10-20 (estimator_modeC_x20.log; table key film_10-25 is the log label). Shared y axis.
(a) three grouped bars per region: energy (1/2 ln k), aligned length L_al, rest (density + strain weight); total as a
    diamond with +-2 SE beside the group, its value printed next to it.
(b) two grouped bars per region: L_eff and alignment ln cos(theta), each with its [min, max] attribution range over
    the 120 orders; L_al = L_eff + align as a diamond (no error bar).
Style: set_mpl_params, Okabe-Ito, no title, no hatching.
"""
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params

DATA = '/home/sk999/Documents/ASME_2027/figures/code/data/'
OUT = '/home/sk999/Documents/ASME_2027/figures/fig_meff_split.pdf'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
REGIONS = [('jet_0-3', 'jet\n' r'$x/D$ 0–3'), ('film_3-10', 'near film\n' r'$x/D$ 3–10'),
           ('film_10-25', 'far film\n' r'$x/D$ 10–20')]
C_EN, C_AL, C_REST, C_LEFF, C_ANG = '#56B4E9', '#E69F00', '#009E73', '#D55E00', '#CC79A7'
PITCH = 1.6


def fmt(v):
    return f'{v:+.2f}'.replace('-', '−')


def read_table(path):
    kind, rows, hdr = '', {}, None
    for ln in open(path, encoding='utf-8'):
        if ln.startswith('# kind:'):
            kind = ln.split(':', 1)[1].strip()
        elif ln.startswith('#') or not ln.strip():
            continue
        elif hdr is None:
            hdr = ln.split()
        else:
            p = ln.split()
            rows[p[0]] = dict(zip(hdr[1:], [float(v) for v in p[1:]]))
    return kind, rows


_, T = read_table(DATA + 'meffR_modeC.txt')
KIND, SE = read_table(DATA + 'meffR_se_time.txt')   # time-batch SE; spatial-block SE kept in meffR_se.txt

print('region | total ± SE | energy L_al rest(rho+Swt) | sum-total | L_eff [min,max] | align [min,max] | L_eff+align-L_al')
for k, _ in REGIONS:
    t = T[k]
    rest = t['density'] + t['Swt']
    print(f"{k:11s} {t['total']:+.3f} ± {SE[k]['SE']:.3f} | {t['energy']:+.3f} {t['L_al']:+.3f} {rest:+.3f} | "
          f"{t['energy'] + t['L_al'] + rest - t['total']:+.4f} | {t['L_eff']:+.3f} [{t['L_eff_min']:+.3f},"
          f"{t['L_eff_max']:+.3f}] | {t['align']:+.3f} [{t['align_min']:+.3f},{t['align_max']:+.3f}] | "
          f"{t['L_eff'] + t['align'] - t['L_al']:+.4f}")

set_mpl_params(0.4, grid_on=True)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.2, 3.3), sharey=True, gridspec_kw=dict(wspace=0.06))

lo = hi = 0.0
# (a)
labs = []
for g, (k, _) in enumerate(REGIONS):
    t, se = T[k], SE[k]['SE']
    x = g * PITCH
    for dx, v, c in ((-0.28, t['energy'], C_EN), (0.0, t['L_al'], C_AL), (0.28, t['density'] + t['Swt'], C_REST)):
        axa.bar(x + dx, v, 0.25, color=c, edgecolor='k', linewidth=0.4)
        lo, hi = min(lo, v), max(hi, v)
    xm, tot = x + 0.62, t['total']
    axa.errorbar(xm, tot, yerr=2 * se, fmt='none', ecolor='k', capsize=1.8, lw=0.8, zorder=5)
    axa.plot(xm, tot, 'D', ms=3.3, mfc='k', mec='k', zorder=6)
    lo, hi = min(lo, tot - 2 * se), max(hi, tot + 2 * se)
    labs.append((xm + 0.07, tot))
# (b)
for g, (k, _) in enumerate(REGIONS):
    t = T[k]
    x = g * PITCH
    for dx, key, c in ((-0.17, 'L_eff', C_LEFF), (0.17, 'align', C_ANG)):
        v, a, b = t[key], t[key + '_min'], t[key + '_max']
        axb.bar(x + dx, v, 0.30, color=c, edgecolor='k', linewidth=0.4)
        axb.errorbar(x + dx, 0.5 * (a + b), yerr=0.5 * (b - a), fmt='none', ecolor='k', capsize=2.0, lw=0.9,
                     zorder=5)
        lo, hi = min(lo, a, v), max(hi, b, v)
    axb.plot(x + 0.55, t['L_al'], 'D', ms=3.3, mfc='k', mec='k', zorder=6)
    lo, hi = min(lo, t['L_al']), max(hi, t['L_al'])

span = hi - lo
axa.set_ylim(lo - 0.08 * span, hi + 0.12 * span)
for xl, tot in labs:
    axa.text(xl, tot, fmt(tot), ha='left', va='center', fontsize=6, zorder=7)

for ax, letter in ((axa, '(a)'), (axb, '(b)')):
    ax.axhline(0, color='k', lw=0.6)
    ax.set_xticks([g * PITCH + 0.15 for g in range(len(REGIONS))])
    ax.set_xticklabels([lab for _, lab in REGIONS], fontsize=6.5)
    ax.set_xlim(-0.6, (len(REGIONS) - 1) * PITCH + 1.05)
    ax.text(0.015, 0.985, letter, transform=ax.transAxes, ha='left', va='top', fontsize=8)
axa.set_ylabel(r'$\ln(\mu_{R,\mathrm{RANS}}/\mu_{R,193})$')

ha = [Patch(facecolor=C_EN, edgecolor='k', linewidth=0.4, label=r'energy $\frac{1}{2}\ln k$'),
      Patch(facecolor=C_AL, edgecolor='k', linewidth=0.4, label=r'aligned length $\ln L_\mathrm{al}$'),
      Patch(facecolor=C_REST, edgecolor='k', linewidth=0.4, label=r'rest ($\ln\rho$ + strain weight)'),
      Line2D([], [], marker='D', color='k', ls='', label=r'total, $\pm2$ SE (' + KIND + ')')]
hb = [Patch(facecolor=C_LEFF, edgecolor='k', linewidth=0.4, label=r'length $\ln L_\mathrm{eff}$'),
      Patch(facecolor=C_ANG, edgecolor='k', linewidth=0.4, label=r'alignment $\ln\cos\theta$'),
      Line2D([], [], marker='D', color='k', ls='', label=r'$\ln L_\mathrm{al}$ (sum)'),
      Line2D([], [], color='k', marker='_', ls='', label='attribution range')]
axa.legend(handles=ha, loc='lower center', ncol=2, bbox_to_anchor=(0.5, 1.0), fontsize=6)
axb.legend(handles=hb, loc='lower center', ncol=2, bbox_to_anchor=(0.5, 1.0), fontsize=6)
fig.savefig(OUT, bbox_inches='tight')
print('wrote', OUT)
