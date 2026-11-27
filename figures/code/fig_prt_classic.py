"""F14 fig_prt_classic: classic turbulent Prandtl number of the reference.

Shows Pr_t = (u''w'' dT/dz) / (w''T'' dU/dz) = eps_m / eps_T (wall-normal form of Kohli and Bogard 2005, Favre
means) of the 193M 3DNS (C193b, stat3 batch b010, 10 FT cumulative mean):
  (a) centre plane y = 0, x/D -2..20, z/D 0..3;
  (b) near-wall strip of (a), x/D 0..20, z/D 0..0.3, vertical axis stretched;
  (c) cross-plane x/D = 10, wall-normal form, y/D -3..3, z/D 0..3;
  (d) cross-plane x/D = 10, lateral form (u''v'' dT/dy) / (v''T'' dU/dy);
  (e) least-squares Pr_t,R(z) = [sum(-rho u''w'' dU/dz)/sum(dU/dz)^2] / [sum(-rho w''T'' dT/dz)/sum(dT/dz)^2] over
      x/D 3-10 (near film) and 10-20 (far film) of the centre plane at each z, z/D 0-0.5; b010 solid with its 95 %
      block-bootstrap band (1 D wide x bins), b006 dashed. Lines: RANS 0.9; Kohli and Bogard (2005): "generally values
      of Pr_t of about 2" at low freestream turbulence (TL 0.6 %), some regions about 0.5; "relatively constant ...
      Pr_t of about 1" at high freestream turbulence (TL 20 %).
Maps: no mask, every point is drawn; white only where the ratio is undefined. Colour bar inferno 0 to 2,
extend='both': negative values take the lowest colour, values above 2 the highest; 0.9 (RANS) is marked. Dashed line
(#56B4E9): inside it both gradients exceed 0.05 U_inf/D (velocity) and 0.05 dT0/D (temperature), dT0 = 100.98 K, as in
the earlier Pr_t,eff map, i.e. where the ratio is well defined.
Printed: median and IQR per region (all points, above thresholds), z range with Pr_t,R >= 1.5, b006 vs b010
agreement (median |ln ratio| per region and for z/D < 0.3, over points where both batches give a positive ratio) and
the freestream turbulence intensity Tu = sqrt(2k/3)/U_inf at x/D -5 and 0, z/D 2-3.
Data: /mnt/Elements/DNS/film_cooling/C193b_BR15, stat3_b010 and stat3_b006, read by prt_corr_common.py (cache
/mnt/Elements/DNS/analysis/asme2027_figs/prt_classes/corr_<full|x10|tu>_b0xx.npz).
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, '/home/sk999/Documents/ASME_2027/figures/code')
from functions import set_mpl_params
import meff_planes_common as mc
import prt_corr_common as cc

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_prt_classic.png'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
VMAX = 2.0
D = mc.D
NBOOT = 1000

set_mpl_params(0.5, grid_on=False)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)

G = dict(np.load(f'{mc.CACHE}/ref193_b010.npz'))


def load(tag, batch, shape):
    f = os.path.join(cc.CORR_DIR, f'corr_{tag}_b{batch:03d}.npz')
    if not os.path.exists(f):
        return None
    z = dict(np.load(f))
    out = {}
    for q in cc.FIELDS:
        a = np.full(int(np.prod(shape)), np.nan, np.float32)
        a[z['idx']] = z[q]
        out[q] = a.reshape(shape).astype(float)
    return out


def classic(c, lateral=False):
    """(Pr_t, well-defined mask) of the wall-normal or the lateral form."""
    if lateral:
        num, dT, den, dU = c['uv'] * c['dTdy'], c['dTdy'], c['vT'] * c['dUdy'], c['dUdy']
    else:
        num, dT, den, dU = c['uw'] * c['dTdz'], c['dTdz'], c['wT'] * c['dUdz'], c['dUdz']
    with np.errstate(divide='ignore', invalid='ignore'):
        v = num / den
    v = np.where(np.isfinite(v), v, np.nan)
    k = (np.abs(dU) > 0.05 * mc.U_INF / D) & (np.abs(dT) > 0.05 * mc.DT0 / D)
    return v, k


def ls_terms(c):
    r = c['rho']
    return np.array([-r * c['uw'] * c['dUdz'], c['dUdz'] ** 2, -r * c['wT'] * c['dTdz'], c['dTdz'] ** 2])


C1 = {b: load('full', b, G['y0_A'].shape) for b in (10, 6)}
SHOW_EXTRA = False     # True: also draw the removed panels (near-wall strip, x/D = 10 cross-planes) and the legend
C2 = {b: (load('x10', b, G['x10_A'].shape) if SHOW_EXTRA else None) for b in (10, 6)}
P1, K1 = classic(C1[10])
have2 = C2[10] is not None
if have2:
    P2, K2 = classic(C2[10])
    P3, K3 = classic(C2[10], lateral=True)
X1, Z1 = G['y0_A'] / D, G['y0_B'] / D
Y2, Z2 = G['x10_A'] / D, G['x10_B'] / D
nz = Z1.shape[0]
zrow = Z1[:, 0]
xrow = X1[0]

W, H = 7.2, (6.2 if SHOW_EXTRA else 4.85)
fig = plt.figure(figsize=(W, H))
kw = dict(cmap='inferno', vmin=0.0, vmax=VMAX, shading='auto', rasterized=True)
ckw = dict(levels=[0.5], colors='#56B4E9', linewidths=0.9, linestyles='--')
LBL = dict(color='black', va='top', bbox=dict(fc='white', ec='none', pad=1.5))
if SHOW_EXTRA:
    ax1 = fig.add_axes([0.7 / W, 5.05 / H, 6.3 / W, 6.3 * 3 / 22 / H])
    axs = fig.add_axes([0.7 / W, 3.75 / H, 6.3 / W, 0.9 / H])
    ax2 = fig.add_axes([0.7 / W, 1.55 / H, 2.0 / W, 1.0 / H])
    ax3 = fig.add_axes([3.0 / W, 1.55 / H, 2.0 / W, 1.0 / H])
    axe = fig.add_axes([5.65 / W, 1.2 / H, 1.35 / W, 2.1 / H])
    cax = fig.add_axes([1.6 / W, 0.40 / H, 2.8 / W, 0.14 / H])
    EX = [(axs, r'(b) $y/D=0$, near wall (stretched)'), (ax2, r'(c) $x/D=10$'), (ax3, r'(d) $x/D=10$, lateral')]
    LAB_A, LAB_E = r'(a) $y/D=0$', '(e)'
else:
    ax1 = fig.add_axes([0.7 / W, 3.85 / H, 6.3 / W, 6.3 * 3 / 22 / H])
    axe = fig.add_axes([0.7 / W, 0.6 / H, 6.3 / W, 1.65 / H])
    cax = fig.add_axes([0.7 / W, 3.1 / H, 6.3 / W, 0.13 / H])
    EX = []
    LAB_A, LAB_E = r'(a) $y/D=0$', '(b)'
m = ax1.pcolormesh(X1, Z1, np.ma.masked_invalid(P1), **kw)
ax1.contour(X1, Z1, K1.astype(float), **ckw)
if SHOW_EXTRA:
    axs.pcolormesh(X1, Z1, np.ma.masked_invalid(P1), **kw)
    axs.contour(X1, Z1, K1.astype(float), **ckw)
    if have2:
        for ax, P, K in ((ax2, P2, K2), (ax3, P3, K3)):
            ax.pcolormesh(Y2, Z2, np.ma.masked_invalid(P), **kw)
            ax.contour(Y2, Z2, K.astype(float), **ckw)
    for ax in (axs, ax2, ax3):
        ax.set_facecolor('white')
        ax.set_ylabel(r'$z/D$')
    for ax in (ax2, ax3):
        ax.set_aspect('equal')
        ax.set_yticks([0, 1, 2, 3])
        ax.set_ylim(0, 3)
    axs.set_xlim(0, 20)
    axs.set_ylim(0, 0.3)
    axs.set_yticks([0, 0.1, 0.2, 0.3])
    ax2.set_xlim(3, -3)                  # y/D increases to the left: plane seen looking downstream
    ax3.set_xlim(3, -3)
    ax2.set_xticks([3, 0, -3])
    ax3.set_xticks([3, 0, -3])
    axs.set_xlabel(r'$x/D$')
    ax2.set_xlabel(r'$y/D$')
    ax3.set_xlabel(r'$y/D$')
    for ax, t in EX:
        ax.text(0.012 if ax is axs else 0.03, 0.95, t, transform=ax.transAxes, **LBL)
ax1.set_facecolor('white')
ax1.set_ylabel(r'$z/D$')
ax1.set_aspect('equal')
ax1.set_yticks([0, 1, 2, 3])
ax1.set_ylim(0, 3)
ax1.set_xlim(-2, 20)
ax1.set_xlabel(r'$x/D$')
ax1.text(0.012, 0.88, LAB_A, transform=ax1.transAxes, **LBL)
if not SHOW_EXTRA:
    ax1.annotate('boundary of the well-defined region', xy=(14.0, 2.07), xytext=(14.0, 1.15),
                 fontsize=8, ha='center', va='center', color='white',
                 arrowprops=dict(arrowstyle='->', color='white', lw=0.8, shrinkA=1, shrinkB=1))
cb = fig.colorbar(m, cax=cax, orientation='horizontal', extend='both')
cb.set_ticks([0, 0.5, 0.9, 1.5, 2.0])
cb.set_ticklabels(['0', '0.5', '0.9', '1.5', '2'])
cb.ax.axvline(0.9, color='white', lw=1.4)
cb.set_label(r'$\mathrm{Pr}_t$')

# (e) x-averaged least-squares profile against z/D, b006 and b010
rng = np.random.default_rng(2)
band = {}
for rn, x0, x1, c, ls in (('near film', 3, 10, '#0072B2', '-'), ('far film', 10, 20, '#009E73', '-')):
    sel = (xrow >= x0) & (xrow < x1)
    binx = np.floor(xrow[sel] - x0).astype(int)
    nb = binx.max() + 1
    for b in (10, 6):
        if C1[b] is None:
            continue
        T = ls_terms(C1[b])[:, :, sel]                         # (4, nz, nx)
        T = np.where(np.isfinite(T), T, 0.0)
        Sb = np.stack([np.stack([T[k][:, binx == j].sum(axis=1) for j in range(nb)], 1) for k in range(4)])  # (4, nz, nb)
        tot = Sb.sum(axis=2)
        with np.errstate(divide='ignore', invalid='ignore'):
            pr = (tot[0] / tot[1]) / (tot[2] / tot[3])
        band[(rn, b)] = pr
        zz = zrow <= 0.5
        if b == 10:
            w = rng.multinomial(nb, np.ones(nb) / nb, size=NBOOT).astype(float)
            with np.errstate(divide='ignore', invalid='ignore'):
                bt = np.stack([(np.einsum('nb,zb->nz', w, Sb[0]) / np.einsum('nb,zb->nz', w, Sb[1])) /
                               (np.einsum('nb,zb->nz', w, Sb[2]) / np.einsum('nb,zb->nz', w, Sb[3]))])[0]
            lo, hi = np.nanpercentile(bt, [2.5, 97.5], axis=0)
            axe.fill_betweenx(zrow[zz], np.clip(lo[zz], -0.5, 3), np.clip(hi[zz], -0.5, 3), color=c, alpha=0.2, lw=0)
            axe.plot(np.clip(pr[zz], -0.5, 3), zrow[zz], color=c, ls='-', label=f'{rn}, b010')
axe.axvline(0.9, color='k', ls=':')
axe.axvline(2.0, color='0.5', ls='-.')
axe.axvline(1.0, color='0.5', ls='--')
axe.set_xlim(-0.5, 3)
axe.set_ylim(0, 0.5)
axe.set_xlabel(r'$\mathrm{Pr}_{t,R}$')
axe.set_ylabel(r'$z/D$')
if SHOW_EXTRA:
    axe.text(0.97, 0.97, LAB_E, transform=axe.transAxes, ha='right', va='top')
else:
    axe.text(0.015, 0.97, LAB_E, transform=axe.transAxes, ha='left', va='top')
if SHOW_EXTRA:
    hh, ll = axe.get_legend_handles_labels()
    fig.legend(hh, ll, loc='lower right', frameon=False, ncol=2, bbox_to_anchor=(1.0, 0.0))
else:
    GR = dict(fontsize=8, ha='left', va='center', arrowprops=dict(arrowstyle='->', lw=0.7, shrinkA=1, shrinkB=1))
    GR['annotation_clip'] = False
    axe.annotate('RANS, 0.9', xy=(0.9, 0.5), xytext=(0.2, 0.56), color='k', **dict(GR, ha='center'))
    axe.annotate('K&B high FST, 1', xy=(1.0, 0.30), xytext=(1.2, 0.30), color='0.35', **GR)
    axe.annotate('K&B low FST, $\\approx2$', xy=(2.0, 0.22), xytext=(2.2, 0.22), color='0.35', **GR)
    hh, ll = axe.get_legend_handles_labels()       # cases only
    axe.legend(hh, ll, loc='upper right', frameon=False, ncol=2, columnspacing=1.2)
fig.savefig(OUT, dpi=300)
print('saved', OUT)


def stats(lab, P, K, s):
    for nm, sel in (('all points', s), ('above thresholds', s & K)):
        v = P[sel]
        v = v[np.isfinite(v)]
        print(f'{lab} [{nm}]: median {np.median(v):.2f}, IQR {np.percentile(v, 25):.2f}-{np.percentile(v, 75):.2f}, '
              f'n {v.size}/{int(s.sum())}')


REG = [('jet x/D 0-3', 0, 3), ('near film 3-10', 3, 10), ('far film 10-20', 10, 20)]
for lab, lo, hi in REG:
    stats(f'centre plane {lab}', P1, K1, (X1 >= lo) & (X1 < hi))
if have2:
    stats('cross-plane x/D=10 wall-normal', P2, K2, np.ones(P2.shape, bool))
    stats('cross-plane x/D=10 lateral', P3, K3, np.ones(P3.shape, bool))
print(f'out of range 0-2 (centre): <0 {np.mean(P1[np.isfinite(P1)] < 0):.3f}, >2 {np.mean(P1[np.isfinite(P1)] > 2):.3f}')

print('\nx-averaged least-squares Pr_t,R(z), z/D 0-0.5')
for (rn, b), pr in sorted(band.items()):
    ok = np.isfinite(pr) & (zrow <= 0.5)
    hi = ok & (pr >= 1.5)
    zs = zrow[hi]
    print(f'{rn} b{b:03d}: z/D with Pr_t,R >= 1.5: ' + (', '.join(f'{z:.2f}' for z in zs) if len(zs) else 'none') +
          f'; at z/D 0.1, 0.2, 0.3, 0.5: ' + ', '.join(f'{pr[np.argmin(abs(zrow - z))]:.2f}' for z in (0.1, 0.2, 0.3, 0.5)))
for rn in ('near film', 'far film'):
    if (rn, 6) in band and (rn, 10) in band:
        a, b = band[(rn, 6)], band[(rn, 10)]
        ok = (zrow <= 0.3) & np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
        print(f'{rn}: median |ln b010/b006| over z/D <= 0.3 rows: {np.median(np.abs(np.log(b[ok] / a[ok]))):.3f} (n rows {ok.sum()})')

# convergence b006 vs b010, pointwise, where both are positive
for tag, Cd, xz, forms in (('centre', C1, (X1, Z1), (False,)), ('cross', C2, (Y2, Z2), (False, True))):
    if Cd[6] is None or Cd[10] is None:
        print(f'{tag}: b006 cache missing, pointwise convergence not computed')
        continue
    for lat in forms:
        a, ka = classic(Cd[6], lat)
        b, kb = classic(Cd[10], lat)
        both = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
        ln = np.abs(np.log(b / np.where(both, a, np.nan)))
        if tag == 'centre':
            sels = [(lab, (xz[0] >= lo) & (xz[0] < hi)) for lab, lo, hi in REG]
            sels += [(lab + ', z/D<0.3', s & (xz[1] < 0.3)) for lab, s in list(sels)]
        else:
            sels = [('x/D=10', np.ones(a.shape, bool)), ('x/D=10, z/D<0.3', xz[1] < 0.3)]
        for lab, s in sels:
            for nm, extra in (('all', True), ('above thresholds', ka & kb)):
                q = both & s & extra
                if q.sum():
                    print(f'convergence {tag} {"lateral" if lat else "wall-normal"} {lab} [{nm}]: median |ln b010/b006| '
                          f'{np.median(ln[q]):.3f}, n {int(q.sum())} (share with both positive {q.sum() / max((s & extra).sum(), 1):.2f})')

# freestream turbulence intensity of the reference
print('\nfreestream turbulence intensity Tu = sqrt(2k/3)/U_inf, z/D 2-3 (b010)')
zt = np.load(os.path.join(cc.CORR_DIR, 'corr_tu_b010.npz'))
v = zt['k'][np.isfinite(zt['k'])].astype(float)
tu = np.sqrt(2 * v / 3) / mc.U_INF
print(f'x/D=-5: median {100 * np.median(tu):.2f} %, range {100 * tu.min():.2f}-{100 * tu.max():.2f} % (n {v.size})')
kk = C1[10]['k']
sel = (np.abs(X1) <= 0.1) & (Z1 >= 2) & (Z1 <= 3) & np.isfinite(kk)
tu0 = np.sqrt(2 * kk[sel] / 3) / mc.U_INF
print(f'x/D=0: median {100 * np.median(tu0):.2f} %, range {100 * tu0.min():.2f}-{100 * tu0.max():.2f} % (n {tu0.size})')
