"""fig_prt_x: classic least-squares Pr_t of the 193M reference (batch b010) as a continuous function of x/D.

Pr_t(x) = nu_R / alpha_R with nu_R = -A/B and alpha_R = -C/E from the span-summed least-squares sums (A, B, C, E of
fig_prt_classic_span.py; 0.02 D bins, summed over the full pitch) of
/mnt/Elements/DNS/analysis/asme2027_figs/n8_lsq/n8_lsq_c193b_b010_jk.npz.
Slabs of x: width 1 D (half-open, [x_c - 0.5, x_c + 0.5)), centres every 0.25 D from 3 to 20; the slab at x/D 20 is
cut by the end of the data (x/D 20), the slab at 3 reaches back to 2.5. Sums are accumulated over the slab and over
the height range of the layer (height mask applied per 0.02 D x bin), then the ratio is formed.
Layers:
  (1) near-wall layer: first bin above the wall (z/D 0.02) up to z_x(x), wall bin z = 0 excluded;
  (2) layer above it: z > z_x(x) up to z/D 2;
  (3) z/D >= 0.3 to the top of the data (z/D 3), for continuity with the regional values of fig_prt_regions.py.
z_x(x): figures/code/data/near_wall_boundary_zx.txt, interpolated linearly onto the x bins (read at run time).
Band: +-2 SE of the delete-one jackknife over the 10 flow-through batches, Pr_t formed per delete-one sample
(jackknife of the ratio), SE = sqrt((n-1)/n sum (theta_b - mean)^2).
Gap: where |alpha_R| <= 2 SE(alpha_R) (jackknife) or Pr_t is not finite, neither line nor band is drawn; the x ranges
are printed. Table: figures/code/data/prt_x_table.txt.
"""
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_prt_x.pdf'
TABLE = '/home/sk999/Documents/ASME_2027/figures/code/data/prt_x_table.txt'
SUMS = '/mnt/Elements/DNS/analysis/asme2027_figs/n8_lsq/n8_lsq_c193b_b010_jk.npz'
ZX = '/home/sk999/Documents/ASME_2027/figures/code/data/near_wall_boundary_zx.txt'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
EPS = 1e-9
W = 1.0                                  # slab width in D
XC = np.arange(3.0, 20.0 + 1e-9, 0.25)   # slab centres
XT = (4, 5, 7, 10, 15, 20)               # table positions

S = np.load(SUMS)
xs, zs = S['x_D'], S['z_D']
Q = {k: S[k] for k in 'ABCEV'}
J = {k: S['jk_' + k] for k in 'ABCEV'}
nb = J['A'].shape[0]

zx_tab = np.loadtxt(ZX)
zx = np.interp(xs, zx_tab[:, 0], zx_tab[:, 1])          # z_x(x) per x bin

ZI = zs[None, :]
LAYERS = [
    ('near-wall layer', '#0072B2', '-', (ZI > 0.01) & (ZI <= zx[:, None] + EPS)),
    ('layer above', '#D55E00', '--', (ZI > zx[:, None] + EPS) & (ZI <= 2.0 + EPS)),
    (r'$z/D\geq0.3$', '#009E73', '-.', np.broadcast_to(ZI >= 0.3 - EPS, (len(xs), len(zs)))),
]


def sums(xmask, hmask, arrs=Q, batch=None):
    m = xmask[:, None] & hmask
    if batch is None:
        return [arrs[k][m].sum() for k in 'ABCE']
    return [arrs[k][batch][m].sum() for k in 'ABCE']


def evaluate(xmask, hmask):
    """Pr_t, SE, alpha_R, SE(alpha_R) for one x mask and one height mask."""
    t = sums(xmask, hmask)
    m = xmask[:, None] & hmask
    jb = np.array([[J[k][b][m].sum() for k in 'ABCE'] for b in range(nb)])
    with np.errstate(divide='ignore', invalid='ignore'):
        pr = (t[0] / t[1]) / (t[2] / t[3])
        al = -t[2] / t[3]
        th = (jb[:, 0] / jb[:, 1]) / (jb[:, 2] / jb[:, 3])
        alb = -jb[:, 2] / jb[:, 3]
    se = np.sqrt((nb - 1) / nb * ((th - th.mean()) ** 2).sum())
    alse = np.sqrt((nb - 1) / nb * ((alb - alb.mean()) ** 2).sum())
    empty = 1.0 - (Q['V'][m] > 0).mean() if m.any() else np.nan
    return pr, se, al, alse, empty


R = []   # per layer: arrays over XC of pr, se, al, alse, empty
for name, col, ls, hm in LAYERS:
    rows = []
    for xc in XC:
        xm = (xs >= xc - W / 2 - EPS) & (xs < xc + W / 2 - EPS) | ((xs >= xc - W / 2 - EPS) & (xc == XC[-1]))
        rows.append(evaluate(xm, hm))
    R.append(np.array(rows))

set_mpl_params(0.5, grid_on=False)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
fig, ax = plt.subplots(figsize=(5.0, 2.8))

# reference lines
for v, lab, ls_, c in ((0.9, 'RANS, 0.9', ':', '#000000'), (0.5, 'free shear, 0.5', '-.', '#CC79A7')):
    ax.axhline(v, ls=ls_, color=c, lw=0.8, zorder=1)
    ax.text(20.25, v, lab, color=c, fontsize=7, ha='left', va='center', clip_on=False)

gaps = {}
for (name, col, ls, hm), r in zip(LAYERS, R):
    pr, se, al, alse = r[:, 0], r[:, 1], r[:, 2], r[:, 3]
    ok = np.isfinite(pr) & np.isfinite(se) & (np.abs(al) > 2 * alse)
    gaps[name] = ~ok
    y = np.where(ok, pr, np.nan)
    ax.fill_between(XC, np.where(ok, pr - 2 * se, np.nan), np.where(ok, pr + 2 * se, np.nan), color=col, alpha=0.2,
                    lw=0, zorder=2)
    ax.plot(XC, y, color=col, ls=ls, zorder=3)

ax.set_xlim(3, 20)
ax.set_xlabel(r'$x/D$')
ax.set_ylabel(r'$\mathrm{Pr}_t$')
ax.set_ylim(0, 1.2)
fig.subplots_adjust(left=0.12, right=0.80, bottom=0.17, top=0.97)
# line labels: text position and arrow tip (x/D, Pr_t)
LABELS = {'near-wall layer': ((7.0, 0.17), (7.0, 0.065)),
          'layer above': ((6.0, 1.05), (6.0, 0.55)),
          r'$z/D\geq0.3$': ((14.0, 0.22), (14.0, 0.43))}
for (name, col, ls, hm) in LAYERS:
    (xt, yt), (xa, ya) = LABELS[name]
    ax.annotate(name, xy=(xa, ya), xytext=(xt, yt), color=col, fontsize=7, ha='center', va='center',
                arrowprops=dict(arrowstyle='->', color=col, lw=0.7, shrinkA=1, shrinkB=1))
fig.savefig(OUT)
print('saved', OUT)

# ---- printed results and table
out = []
out.append('# Classic least-squares Pr_t(x) of the 193M reference b010: slabs 1 D wide, centres every 0.25 D, span-summed.')
out.append('# Pr_t = nu_R/alpha_R from summed numerators/denominators; SE = delete-one jackknife (10 batches) of the ratio.')
out.append('# layer 1: z/D 0.02..z_x(x); layer 2: z_x(x)..2; layer 3: z/D >= 0.3 to 3. empty = share of empty 0.02 D bins in slab x layer.')
out.append('# x/D  layer  Pr_t  SE  alpha_R  SE(alpha_R)  empty  resolved')
for xt in XT:
    i = int(np.argmin(np.abs(XC - xt)))
    for li, ((name, *_), r) in enumerate(zip(LAYERS, R), 1):
        p, s, a, as_, e = r[i]
        out.append(f'{xt:5.1f}  {li}  {p:7.4f}  {s:7.4f}  {a:10.3e}  {as_:10.3e}  {e:5.2f}  {"yes" if abs(a) > 2 * as_ else "no"}')
out.append('# gaps (alpha_R within 2 SE of zero or Pr_t not finite), x/D ranges of slab centres:')
for name, g in gaps.items():
    rng_ = []
    i = 0
    while i < len(g):
        if g[i]:
            j = i
            while j + 1 < len(g) and g[j + 1]:
                j += 1
            rng_.append(f'{XC[i]:.2f}-{XC[j]:.2f}')
            i = j + 1
        else:
            i += 1
    out.append(f'#   {name}: ' + (', '.join(rng_) if rng_ else 'none'))
out.append('# layer 3 regional fits (one fit over the region, same convention as fig_prt_regions.py, z/D >= 0.3):')
for lo, hi, ref in ((3, 10, 0.225), (10, 20, 0.454)):
    xm = (xs >= lo - EPS) & (xs < hi - EPS)
    hm = LAYERS[2][3]
    p, s, *_ = evaluate(xm, hm)
    out.append(f'#   x/D {lo}-{hi}: Pr_t {p:.4f}  SE {s:.4f}  (expected {ref})')
    sel = (XC >= lo) & (XC <= hi)
    out.append(f'#   mean of the slab values over centres {lo}-{hi}: {np.nanmean(R[2][sel, 0]):.4f} (for information)')
out.append('# empty-bin fraction per slab (layers 1, 2, 3) at the table positions is in the column "empty".')
txt = '\n'.join(out)
print(txt)
with open(TABLE, 'w') as f:
    f.write(txt + '\n')
