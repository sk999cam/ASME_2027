"""F14 line plot: exact factor split of the eddy-viscosity difference in the near-wall layer against x/D.

Usage: python3 fig_meff_split_film.py     (needs data/near_wall_boundary_zx.txt, written by fig_meff_maps.py)
Output: /home/sk999/Documents/ASME_2027/figures/fig_meff_split_film.pdf

Near-wall layer = wall to z_x(x), the height where the ratio of the conditional span-mean eddy viscosity RANS/193M first reaches 1 going up
from the wall (curve read from data/near_wall_boundary_zx.txt). Wall nodes excluded, z/D <= 2. Cells with the point-wise 193M |S*| D/U_inf > 0.05
are kept (3D mask), everything inside is a plain volume statistic. Evaluated in slabs of x of width 1 D (|x - centre| < 0.5 D), centres every
0.25 D from 4.5 to 19.5 (x/D 4-20). z_x from the full-pitch conditional span means (meff_span files); the split itself is from the point file
(exact products per point). A bin-mean cross-check from the meff_span files is printed.
(a) total = ln( sum V mu_t,RANS / sum V mu_eff,193M ) and its exact Shapley split into rho, k^(1/2), L_eff, cos(theta) (four players).
(b) the ratios RANS/193M of the volume-mean k^(1/2), L (RANS L_R = mu_t/(rho k^(1/2)), 193M L_eff) and 1/cos(theta) (RANS cos = 1).
The split is computed in _meff_split_nearfar.compute_curve.
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from matplotlib.lines import Line2D

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from functions import set_mpl_params
import _meff_split_nearfar as SP

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_meff_split_film.pdf'
PNG = '/home/sk999/Documents/ASME_2027/planning/review/56_fig_meff_split_les.png'    # review render (200 dpi)
ZXF = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'near_wall_boundary_zx.txt')
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']

zt = np.loadtxt(ZXF)
CEN = np.arange(4.5, 19.5 + 1e-9, 0.25)         # slabs of 1 D: x/D 4-20 (first slab 4-5; owner 2026-11-27: z_x confirmed from x/D 4)
FMIN = float(os.environ.get('F14_FMIN', 0.2))
ZXT = (zt[:, 0], zt[:, 1])
# primary: point file (pts.npz), 3D strain mask, products exact per point; layer = wall to the full-pitch z_x
cc_ = lambda c_, w_=1.0, **kw: SP.compute_curve(ZXT, c_, w_, **kw)  # noqa: E731
cur = cc_(CEN)
x = np.array([r['x'] for r in cur])
tot = np.array([r['tot'] for r in cur])
p4 = np.array([r['p4'] for r in cur])
res = max(abs(r['res']) for r in cur)
print(f'point-file split, {len(cur)} slabs; max |sum of the four terms - total| = {res:.3e}')
assert res < 0.005
fm = cc_([12.0], 16.0)[0]
print(f'first slab (x/D 4-5, centre {cur[0]["x"]}): total {cur[0]["tot"]:+.3f}; rho, k^1/2, L_eff, cos: ' + ' '.join(f'{v:+.3f}' for v in cur[0]['p4']) +
      f'; kept {cur[0]["kept"]:.3f}')
print(f'near-wall layer mean x/D 4-20: total {fm["tot"]:+.3f}; rho, k^1/2, L_eff, cos: ' + ' '.join(f'{v:+.3f}' for v in fm['p4']) +
      f'; gap {fm["res"]:+.1e}')
fa = cc_([12.0], 16.0, layer='above')[0]
va = fa['vmean']
print(f'above z_x, x/D 4-20 (point file): mu ratio {np.exp(fa["tot"]):.3f} k^1/2 {va["k"][0] / va["k"][1]:.3f} L {va["L"][0] / va["L"][1]:.3f} '
      f'1/cos {1 / va["c"][1]:.3f}')
cura = cc_(CEN, layer='above')
for nm_, fn_ in (('mu', lambda r: np.exp(r['tot'])), ('k^1/2', lambda r: r['vmean']['k'][0] / r['vmean']['k'][1]),
                 ('L_eff', lambda r: r['vmean']['L'][0] / r['vmean']['L'][1]), ('1/cos', lambda r: 1 / r['vmean']['c'][1])):
    v_ = np.array([fn_(r) for r in cura])
    print(f'  above z_x range over slabs x/D 4-20, {nm_}: {v_.min():.3f}-{v_.max():.3f}')
for xq in (5, 10, 15, 20):
    r = cc_([xq], 1.0, layer='above')[0]
    vm = r['vmean']
    print(f'  above z_x, slab x/D {xq}: mu ratio {np.exp(r["tot"]):.3f} k^1/2 {vm["k"][0] / vm["k"][1]:.3f} L {vm["L"][0] / vm["L"][1]:.3f} 1/cos {1 / vm["c"][1]:.3f}')
# cross-check: bin means of the meff_span files (gap = mean of a product against product of the means)
import _meff_span_split as SS2
SG = SS2.load()
xcur = SS2.compute_curve(SG, ZXT, CEN, 1.0, fmin=FMIN)
print(f'cross-check, bin means (f_min {FMIN}): max |gap| {max(abs(r["res"]) for r in xcur):.3f}; totals at x/D 5, 10, 19.5: ' +
      ', '.join(f'{xcur[int(np.abs(x - q).argmin())]["tot"]:+.3f}' for q in (5, 10, 19.5)))
for xq in (5, 10, 20):
    j = int(np.abs(x - min(xq, 19.5)).argmin())
    print(f'x/D {x[j]:.2f}: total {tot[j]:+.3f}; rho, k^1/2, L_eff, cos: ' + ' '.join(f'{v:+.3f}' for v in p4[j]))
# 20M LES (resolution ladder, owner 2026-11-27): total only, from the meff_span files (the point file has no LES). Route: conditional
# Delta-y-weighted span means cy_mu_mean (span_C20bW.npz against span_C193b.npz), rebinned to 0.1 D in x and put on the 193M z levels
# (_meff_span_split.load), weights w = dx dz f (f = 193M valid fraction >= f_min), layer wall to z_x (and above z_x), slabs of 1 D.
# RANS by the same route is computed for comparison. Table: data/meff_les_table.txt.
LES = os.environ.get('F14_LES', '1') == '1'
CY = dict(mu='cy_mu_mean', rho='cy_rho_mean', k='cy_sqrtk_mean', L='cy_Leff_mean', c='cy_cos_mean')
if LES:
    _kc = dict(SS2.KEYS)
    SS2.KEYS = CY
    try:
        GL, GRc = SS2.load(case='C20bW'), SS2.load(case='RANS50')
    finally:
        SS2.KEYS = _kc
    les_n, les_a = SS2.compute_curve(GL, ZXT, CEN, 1.0, fmin=FMIN), SS2.compute_curve(GL, ZXT, CEN, 1.0, fmin=FMIN, layer='above')
    rsp_n, rsp_a = SS2.compute_curve(GRc, ZXT, CEN, 1.0, fmin=FMIN), SS2.compute_curve(GRc, ZXT, CEN, 1.0, fmin=FMIN, layer='above')
    les_x = np.array([r['x'] for r in les_n])
    assert np.allclose(les_x, x)
    tot_les = np.array([r['tot'] for r in les_n])
    mean_ = lambda G, lay: SS2.compute_curve(G, ZXT, [12.0], 16.0, fmin=FMIN, layer=lay)[0]['tot']  # noqa: E731
    XQT = (5.0, 10.0, 19.5)
    pick = lambda cur_, xq: cur_[int(np.abs(np.array([r['x'] for r in cur_]) - xq).argmin())]['tot']  # noqa: E731
    PAPER_N = {5.0: -0.415, 10.0: -0.548, 19.5: -0.589}
    T = ['# 20M LES against the 193M reference (b010): eddy-viscosity ratio in the near-wall layer (below z_x) and above z_x, x/D 4-20',
         '# ratio = sum(w mu_X) / sum(w mu_193), mu = cy_mu_mean (conditional, Delta-y-weighted span means, 193M |S*| D/U_inf > 0.05 mask point by point),',
         f'# span_<case>.npz of meff_span rebinned to 0.1 D in x, on the 193M z levels; w = dx dz f, f = 193M valid fraction, f >= {FMIN}; z/D <= 2, wall row excluded;',
         '# slabs of 1 D in x (|x - centre| < 0.5 D); layer means over x/D 4-20 (one slab 16 D wide). z_x = data/near_wall_boundary_zx.txt.',
         '# LES mu_eff = resolved-stress effective viscosity + mean WALE stress. "RANS span" = RANS 50M by the same route as the LES (comparison).',
         '# "RANS point file" = RANS 50M by the paper route (_meff_split_nearfar: point-wise 3D mask, mesh-volume weights, 193M interpolated to the RANS points).',
         '', 'BELOW z_x (near-wall layer): ln ratio / ratio',
         f'  {"x/D":>9s} {"20M LES":>16s} {"RANS span":>16s} {"RANS point file":>16s} {"paper":>8s} {"RANS span - paper":>18s}']
    for xq in XQT:
        l_, r_, p_ = pick(les_n, xq), pick(rsp_n, xq), pick(cur, xq)
        T.append(f'  {xq:9.1f} {l_:+8.3f} / {np.exp(l_):5.3f} {r_:+8.3f} / {np.exp(r_):5.3f} {p_:+8.3f} / {np.exp(p_):5.3f} {PAPER_N[xq]:+8.3f} {r_ - PAPER_N[xq]:+18.3f}')
    l_, r_, p_ = mean_(GL, 'near'), mean_(GRc, 'near'), fm['tot']
    T.append(f'  {"mean 4-20":>9s} {l_:+8.3f} / {np.exp(l_):5.3f} {r_:+8.3f} / {np.exp(r_):5.3f} {p_:+8.3f} / {np.exp(p_):5.3f}')
    T += ['', 'ABOVE z_x (to the strain-mask edge, z/D <= 2): ln ratio / ratio',
          f'  {"x/D":>9s} {"20M LES":>16s} {"RANS span":>16s} {"RANS point file":>16s} {"paper (RANS)":>13s}']
    for xq in XQT:
        l_, r_, p_ = pick(les_a, xq), pick(rsp_a, xq), pick(cura, xq)
        T.append(f'  {xq:9.1f} {l_:+8.3f} / {np.exp(l_):5.3f} {r_:+8.3f} / {np.exp(r_):5.3f} {p_:+8.3f} / {np.exp(p_):5.3f} {"1.00-1.12":>13s}')
    l_, r_, p_ = mean_(GL, 'above'), mean_(GRc, 'above'), fa['tot']
    T.append(f'  {"mean 4-20":>9s} {l_:+8.3f} / {np.exp(l_):5.3f} {r_:+8.3f} / {np.exp(r_):5.3f} {p_:+8.3f} / {np.exp(p_):5.3f}')
    la, ra = np.exp([r['tot'] for r in les_a]), np.exp([r['tot'] for r in rsp_a])
    T += ['', f'Range of the ratio over the slabs x/D 4.5-19.5 above z_x: 20M LES {la.min():.3f}-{la.max():.3f}; RANS (same route) {ra.min():.3f}-{ra.max():.3f}',
          f'Range below z_x: 20M LES {np.exp(tot_les).min():.3f}-{np.exp(tot_les).max():.3f}; RANS (same route) {np.exp([r["tot"] for r in rsp_n]).min():.3f}-'
          f'{np.exp([r["tot"] for r in rsp_n]).max():.3f}']
    open(os.path.join(os.path.dirname(ZXF), 'meff_les_table.txt'), 'w').write('\n'.join(T) + '\n')
    print('\n'.join(T))
rk = np.array([r['vmean']['k'][0] / r['vmean']['k'][1] for r in cur])
rl = np.array([r['vmean']['L'][0] / r['vmean']['L'][1] for r in cur])
rc = np.array([1.0 / r['vmean']['c'][1] for r in cur])
# same three ratios for the strained film above z_x (same point-wise method, layer='above'), dashed in panel (b)
rka = np.array([r['vmean']['k'][0] / r['vmean']['k'][1] for r in cura])
rla = np.array([r['vmean']['L'][0] / r['vmean']['L'][1] for r in cura])
rca = np.array([1.0 / r['vmean']['c'][1] for r in cura])
assert np.allclose([r['x'] for r in cura], x)
for nm_, v_ in (('k^1/2', rka), ('L', rla), ('1/cos', rca)):
    m_ = x >= 5
    print(f'  above z_x, panel (b), x/D 5-20, {nm_}: {v_[m_].min():.3f}-{v_[m_].max():.3f}')

names = ['total', r'$\rho$', r'$k^{1/2}$', r'$L_\mathrm{eff}$', r'$\cos\theta$']
cols = ['#000000', '#CC79A7', '#56B4E9', '#009E73', '#D55E00']      # rho pink: #E69F00 dash-dot is the 20M LES (paper-wide)
styles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1, 1, 1))]

# Drawn at the printed size for a full column (3.56 in; include width=\columnwidth, scale 1); body text 10 pt -> sf = 10/22.
set_mpl_params(10 / 22, grid_on=False)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
fig, (ax, ax2, ax3) = plt.subplots(1, 3, figsize=(7.5, 2.55))
series = [tot] + [p4[:, i] for i in range(4)]
# direct labels inside (a): (offset above or below the line end)
voff = {'total': 'bottom', r'$\rho$': 'bottom', r'$k^{1/2}$': 'bottom', r'$L_\mathrm{eff}$': 'top', r'$\cos\theta$': 'bottom'}
if LES:
    voff[r'$\rho$'] = 'top'                 # the LES line runs above the rho line at the right end; the rho label goes below its line
for s, n_, c_, ls in zip(series, names, cols, styles):
    ax.plot(x, s, color=c_, ls=ls, lw=1.5 if n_ == 'total' else 0.8)
    ax.text(19.8, s[-1] + (0.02 if voff[n_] == 'bottom' else -0.02), n_, color=c_, ha='right', va=voff[n_])
if LES:
    ax.plot(x, tot_les, color='#E69F00', ls='-.', lw=1.5, label='20M LES')            # total only, not the factor split
    ax.legend(loc='upper left', frameon=False, handlelength=2.2, borderaxespad=0.3, fontsize=plt.rcParams['font.size'] - 1)
ax.axhline(0, color='k', lw=0.5)
ax.set_xlim(4, 20)
ax.set_ylim(-0.85, 0.75 if LES else 0.6)
ax.set_ylabel(r'$\ln(\mu/\mu_{193\mathrm{M}})$' if LES else r'$\ln(\mu_\mathrm{RANS}/\mu_{193\mathrm{M}})$')
fac = ((cols[2], styles[2], r'$k^{1/2}$'), (cols[3], styles[3], r'$L$'), (cols[4], styles[4], r'$1/\cos\theta$'))
for a_, vals, lab in ((ax2, (rk, rl, rc), r'below $z_x$'), (ax3, (rka, rla, rca), r'above $z_x$')):
    for s, (c_, ls, _) in zip(vals, fac):
        a_.plot(x, s, color=c_, ls=ls, lw=0.8)
    a_.axhline(1, color='k', lw=0.5)
    a_.set_ylim(0, 3.2)
    a_.set_xlim(4, 20)
    a_.text(0.04, 0.96, lab, transform=a_.transAxes, ha='left', va='top')
ax2.set_ylabel('ratio RANS / 193M')
ax2.legend(handles=[Line2D([], [], color=c_, ls=ls, lw=0.8, label=l_) for c_, ls, l_ in fac], loc='upper right',
           fontsize=plt.rcParams['font.size'] - 1, handlelength=2.2, borderaxespad=0.3, labelspacing=0.2, frameon=False)
for a_, t in ((ax, '(a)'), (ax2, '(b)'), (ax3, '(c)')):
    a_.set_xlabel(r'$x/D$')
    a_.set_xticks([5, 10, 15, 20])
    a_.text(0.0, 1.02, t, transform=a_.transAxes, ha='left', va='bottom')
fig.tight_layout(pad=0.3, w_pad=0.6)
print('axes box (w, h) in:', [(round(a_.get_position().width * 7.5, 2), round(a_.get_position().height * 2.55, 2)) for a_ in (ax, ax2, ax3)])
fig.savefig(OUT)
fig.savefig(PNG, dpi=200)
print('saved', OUT)
