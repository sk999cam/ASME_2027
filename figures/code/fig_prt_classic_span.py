"""fig_prt_classic_span: span-averaged version of fig_prt_classic (comparison file, centre-plane figure unchanged).

Pr_t(x, z) = (A/B) / (C/E) with the y-summed least-squares sums of the 193M 3DNS (stat3 b010), file
/mnt/Elements/DNS/analysis/asme2027_figs/n8_lsq/n8_lsq_c193b_b010.npz (0.02 D bins in x and z; written by
4dns_filmcooled/analysis/eps_prt_c193b/n8_sums.py):
    A = sum_y V rho u''w'' dU/dz,  B = sum_y V (dU/dz)^2,  C = sum_y V rho w''T'' dT/dz,  E = sum_y V (dT/dz)^2.
The ratio (A/B)/(C/E) needs no sign change and no density.
  (a) map over x/D -2..20, z/D 0..3. A, B, C, E are accumulated onto a regular grid of 0.1 D in x and 0.04 D in z
      (5 x 2 original bins, no empty cell) and the ratio is formed afterwards; nothing is interpolated or masked.
      Dashed line (#56B4E9): well-defined region, where the rms gradients of the cell, sqrt(B/V) and sqrt(E/V),
      exceed 0.05 U_inf/D (velocity) and 0.05 dT0/D (temperature), dT0 = 100.98 K; the threshold of the centre-plane
      figure applied to the span-rms gradient (V = sum of the grid-point volumes in the cell).
  (b) Pr_t against z/D from the sums over x/D 3-10 (near film), 10-20 (far film) and 0-3 (jet, thin dotted) at 0.02 D
      in z; thin dashed lines of the same colours: the centre-plane b010 profiles of fig_prt_classic.
      Optional 95 % band: pass a per-batch file (npz with A, B, C, E of shape [nbatch, nx, nz] on the grid of the
      sums file) as argv[1] or in BATCH_FILE; the band is the percentile interval of a batch bootstrap.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler
from scipy import ndimage
import matplotlib.ticker  # noqa: F401

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, '/home/sk999/Documents/ASME_2027/figures/code')
from functions import set_mpl_params
import meff_planes_common as mc
import prt_corr_common as cc

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_prt_classic_span.png'
SUMS = '/mnt/Elements/DNS/analysis/asme2027_figs/n8_lsq/n8_lsq_c193b_b010.npz'
BATCH_FILE = sys.argv[1] if len(sys.argv) > 1 else os.environ.get('BATCH_FILE', '')
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
VMAX = 2.0
D = mc.D
NBOOT = 1000
CX, CZ = 5, 2                       # original bins per map cell in x and z (0.1 D x 0.04 D)
ZMAX_B = 1.0                     # upper z/D limit of panel (b)

set_mpl_params(0.5, grid_on=False)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)

S = np.load(SUMS)
xs, zs = S['x_D'], S['z_D']
Q = {k: S[k] for k in ('A', 'B', 'C', 'E', 'V')}


def ratio(A, B, C, E):
    with np.errstate(divide='ignore', invalid='ignore'):
        return (A / B) / (C / E)


# ---- (a) coarse grid x/D -2..20 (bins centred -2.00..19.98), z/D 0..3 (bins centred 0.00..2.98 plus the last pair)
nxb = int(round((20.0 - (-2.0)) / 0.02))              # 1100 original bins in x
nzb = (len(zs) // CZ) * CZ                              # 150 bins in z (z/D 0 to 2.98); the 3.00 bin is not used
xe = xs[:nxb].reshape(-1, CX).mean(1)
ze = zs[:nzb].reshape(-1, CZ).mean(1)


def coarse(M):
    return M[:nxb, :nzb].reshape(nxb // CX, CX, nzb // CZ, CZ).sum(axis=(1, 3))


Cg = {k: coarse(v) for k, v in Q.items()}
P1 = ratio(Cg['A'], Cg['B'], Cg['C'], Cg['E'])
with np.errstate(divide='ignore', invalid='ignore'):
    K1 = (np.sqrt(Cg['B'] / Cg['V']) > 0.05 * mc.U_INF / D) & (np.sqrt(Cg['E'] / Cg['V']) > 0.05 * mc.DT0 / D)
# island removal: connected pieces (4-neighbour) of the well-defined region, and of its complement, smaller than
# MIN_AREA (in D^2) are merged into the surrounding class before the boundary is drawn
MIN_AREA = 0.1
K_RAW = K1.copy()
removed = []
for val in (True, False):
    lab, nl = ndimage.label(K1 == val)
    for i in range(1, nl + 1):
        m_ = lab == i
        area = m_.sum() * 0.1 * 0.04
        if area < MIN_AREA:
            ii, jj = np.nonzero(m_)
            removed.append((val, area, xe[ii].mean(), ze[jj].mean()))
            K1[m_] = not val
# cell edges for pcolormesh
dx, dz = 0.02 * CX, 0.02 * CZ
XE = np.concatenate([xe - dx / 2, [xe[-1] + dx / 2]])
ZE = np.concatenate([ze - dz / 2, [ze[-1] + dz / 2]])
ZE[0] = 0.0
XEg, ZEg = np.meshgrid(XE, ZE, indexing='ij')
XC, ZC = np.meshgrid(xe, ze, indexing='ij')

# ---- wall units: u_tau and nu_w per region from the span- and x-averaged mean wall shear of the 193M b010
WALL = '/mnt/Elements/DNS/analysis/asme2027_figs'
_b = np.load(f'{WALL}/bl_extract/c193b_b010.npz')        # wall columns: x/D, y/D, C_f, u_e/U_inf
_p = np.load(f'{WALL}/fig_eta_wall/plate_193M_b010.npz')  # same points: x/D, y/D, eta
RHO_INF, P_EXIT, R_GAS = 0.88918, 76596.2, 287.143
_Tw = 300.0 - _p['e'] * mc.DT0
_c1 = 2.5e-4 * 296.296 ** -1.5 * (296.296 + 110.4)
_mu = _c1 * _Tw ** 1.5 / (_Tw + 110.4)
_rho = P_EXIT / (R_GAS * _Tw)
_tau = _b['cf'] * 0.5 * RHO_INF * (_b['ue'] * mc.U_INF) ** 2


def wall_units(lo, hi):
    s_ = (_b['x'] >= lo) & (_b['x'] < hi) & (_b['cf'] > 0)
    inv = np.unique(np.floor(_b['x'][s_] * 10).astype(int) * 100 + np.floor((_b['y'][s_] + 3) * 5).astype(int),
                    return_inverse=True)[1]
    cnt = np.bincount(inv)
    mean = lambda a: (np.bincount(inv, a[s_]) / cnt).mean()     # equal weight per 0.1 D x 0.2 D cell
    tau, rho, mu = mean(_tau), mean(_rho), mean(_mu)
    return dict(tau=tau, rho=rho, nu=mu / rho, ut=np.sqrt(tau / rho), cf=mean(_b['cf']))


# ---- centre-plane profiles (existing figure), b010
G = dict(np.load(f'{mc.CACHE}/ref193_b010.npz'))
zrow, xrow = G['y0_B'][:, 0] / D, G['y0_A'][0] / D
f = os.path.join(cc.CORR_DIR, 'corr_full_b010.npz')
zz_ = dict(np.load(f))
c = {}
for q in cc.FIELDS:
    a = np.full(int(np.prod(G['y0_A'].shape)), np.nan, np.float32)
    a[zz_['idx']] = zz_[q]
    c[q] = a.reshape(G['y0_A'].shape).astype(float)
r_ = c['rho']
LT = np.array([-r_ * c['uw'] * c['dUdz'], c['dUdz'] ** 2, -r_ * c['wT'] * c['dTdz'], c['dTdz'] ** 2])
LT = np.where(np.isfinite(LT), LT, 0.0)

REG = [('near film', 3, 10, '#0072B2', '-'), ('far film', 10, 20, '#009E73', '-')]
WU = {rn: wall_units(lo, hi) for rn, lo, hi, _, _ in REG}
LS = {rn: WU[rn]['nu'] / WU[rn]['ut'] / D for rn in WU}          # viscous length / D


def flags(grad_ok, nu, al, z):
    """True = Pr_t not defined: gradient test failed, counter-gradient (nu_R <= 0 or alpha_R <= 0), or the bin on either
    side of a sign change of nu_R or alpha_R (the ratio passes through a pole there)."""
    bad = ~grad_ok | ~(nu > 0) | ~(al > 0)
    for q in (nu, al):
        sg = np.sign(q)
        ch = np.nonzero((sg[1:-1] != sg[2:]) & np.isfinite(q[1:-1]) & np.isfinite(q[2:]))[0] + 1   # z > 0 bins only
        bad[ch] = True
        bad[ch + 1] = True
    return bad


def span_profile(lo, hi):
    sel = (xs >= lo - 1e-9) & (xs < hi - 1e-9)
    t = [Q[k][sel].sum(0) for k in 'ABCE']
    Vs = Q['V'][sel].sum(0)
    with np.errstate(divide='ignore', invalid='ignore'):
        gok = (np.sqrt(t[1] / Vs) > 0.05 * mc.U_INF / D) & (np.sqrt(t[3] / Vs) > 0.05 * mc.DT0 / D)
        nu, al = -t[0] / t[1], -t[2] / t[3]                  # nu_R = -A/B, alpha_R = -C/E (positive in the film)
    return ratio(*t), flags(gok, nu, al, zs), nu, al


def centre_profile(lo, hi):
    sel = (xrow >= lo) & (xrow < hi)
    t = LT[:, :, sel].sum(2)                                  # LT already holds -rho uw dU/dz and -rho wT dT/dz
    with np.errstate(divide='ignore', invalid='ignore'):
        gok = (np.sqrt(t[1] / sel.sum()) > 0.05 * mc.U_INF / D) & (np.sqrt(t[3] / sel.sum()) > 0.05 * mc.DT0 / D)
        nu, al = t[0] / t[1], t[2] / t[3]
    return (t[0] / t[1]) / (t[2] / t[3]), flags(gok, nu, al, zrow), nu, al


# optional per-batch band
band = {}
if BATCH_FILE and os.path.exists(BATCH_FILE):
    bz = np.load(BATCH_FILE)
    rng = np.random.default_rng(2)
    nb = bz['A'].shape[0]
    for rn, lo, hi, _, _ in REG:
        sel = (xs >= lo - 1e-9) & (xs < hi - 1e-9)
        sb = np.stack([bz[k][:, sel].sum(1) for k in 'ABCE'])          # (4, nb, nz)
        w = rng.multinomial(nb, np.ones(nb) / nb, size=NBOOT).astype(float)
        t = [np.einsum('nb,bz->nz', w, sb[i]) for i in range(4)]
        bt = ratio(*t)
        band[rn] = np.nanpercentile(bt, [2.5, 97.5], axis=0)

W, H = 7.2, 5.4
fig = plt.figure(figsize=(W, H))
kw = dict(cmap='inferno', vmin=0.0, vmax=VMAX, shading='auto', rasterized=True)
ckw = dict(levels=[0.5], colors='#56B4E9', linewidths=0.9, linestyles='--')
LBL = dict(color='black', va='top', bbox=dict(fc='white', ec='none', pad=1.5))
ax1 = fig.add_axes([0.7 / W, 4.4 / H, 6.3 / W, 6.3 * 3 / 22 / H])
axe = fig.add_axes([0.7 / W, 0.6 / H, 5.9 / W, 2.2 / H])
cax = fig.add_axes([0.7 / W, 3.65 / H, 6.3 / W, 0.13 / H])

m = ax1.pcolormesh(XEg, ZEg, np.ma.masked_invalid(P1), **kw)
ax1.contour(XC, ZC, K1.astype(float), **ckw)
ax1.set_facecolor('white')
ax1.set_ylabel(r'$z/D$')
ax1.set_aspect('equal')
ax1.set_yticks([0, 1, 2, 3])
ax1.set_ylim(0, 3)
ax1.set_xlim(-2, 20)
ax1.set_xlabel(r'$x/D$')
ax1.text(0.012, 0.88, r'(a) span-averaged', transform=ax1.transAxes, **LBL)
ax1.annotate('boundary of the well-defined region', xy=(14.0, 2.07), xytext=(14.0, 1.15),
             fontsize=8, ha='center', va='center', color='white',
             arrowprops=dict(arrowstyle='->', color='white', lw=0.8, shrinkA=1, shrinkB=1))
cb = fig.colorbar(m, cax=cax, orientation='horizontal', extend='both')
cb.set_ticks([0, 0.5, 0.9, 1.5, 2.0])
cb.set_ticklabels(['0', '0.5', '0.9', '1.5', '2'])
cb.ax.axvline(0.9, color='white', lw=1.4)
cb.set_label(r'$\mathrm{Pr}_t$')

# ---- (b) Pr_t against z+ (log), each region with its own u_tau; z/D axis on the right
prof = {}
zp0 = zs > 0                                   # the wall bin z = 0 has no z+
XL = (-0.5, 3.0)


def draw(axe, p, bad, zp, col, ls, lw, ms, sel):
    """Unflagged bins as the full line, flagged bins faint with open markers; flagged values outside the x range
    are drawn at the axis edge with a small arrow."""
    axe.plot(np.where(bad, np.nan, p)[sel], zp[sel], color=col, ls=ls, lw=lw)
    axe.plot(np.where(bad & (p >= XL[0]) & (p <= XL[1]), p, np.nan)[sel], zp[sel], color=col, ls=ls, lw=0.4, alpha=0.5)
    axe.plot(np.where(bad, p, np.nan)[sel], zp[sel], color=col, ls='none', marker='o', mfc='none', ms=ms, mew=0.5, alpha=0.7)
    for v, z in zip(p[sel & bad], zp[sel & bad]):
        if np.isfinite(v) and (v > XL[1] or v < XL[0]):
            edge, d = (XL[1], -1) if v > XL[1] else (XL[0], 1)
            axe.annotate('', xy=(edge, z), xytext=(edge + d * 0.18, z), textcoords='data',
                         arrowprops=dict(arrowstyle='->', color=col, lw=0.6, alpha=0.7, shrinkA=0, shrinkB=0),
                         annotation_clip=False)


FLAG = {}
for rn, lo, hi, col, ls in REG:
    ps, bad, nu, al = span_profile(lo, hi)
    pc, badc, nuc, alc = centre_profile(lo, hi)
    zp = zs / LS[rn]
    zc = zrow / LS[rn]
    prof[rn] = (ps, bad, pc, zp)
    FLAG[rn] = (bad, badc, nu, al, nuc, alc)
    sel = zp0 & (zs <= ZMAX_B + 1e-9)
    draw(axe, ps, bad, zp, col, '-', 1.0, 2.5, sel)
    zcs = (zrow > 0) & (zrow <= ZMAX_B + 1e-9)
    draw(axe, pc, badc, zc, col, '--', 0.6, 1.8, zcs)
    if rn in band:
        lo_, hi_ = band[rn]
        axe.fill_betweenx(zp[sel], lo_[sel], hi_[sel], color=col, alpha=0.2, lw=0)
axe.set_yscale('log')
axe.axvline(0.9, color='k', ls=':')
axe.axvline(2.0, color='0.5', ls='-.')
axe.axvline(1.0, color='0.5', ls='--')
axe.set_xlim(*XL)
ZP_LO, ZP_HI = 2.0, ZMAX_B / LS['near film']
axe.set_ylim(ZP_LO, ZP_HI)
axe.set_xlabel(r'$\mathrm{Pr}_{t,R}$')
axe.set_ylabel(r'$z^+$')
axe.text(0.015, 0.97, '(b)', transform=axe.transAxes, ha='left', va='top')
sec = axe.secondary_yaxis('right', functions=(lambda zp: zp * LS['near film'], lambda zd: zd / LS['near film']))
sec.set_yscale('log')
sec.set_yticks([0.02, 0.05, 0.1, 0.2, 0.5, 1.0])
sec.set_yticklabels(['0.02', '0.05', '0.1', '0.2', '0.5', '1'])
sec.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
sec.set_ylabel(r'$z/D$ (near-film $u_\tau$)')
axe.set_yticks([2, 5, 10, 20, 50, 100])
axe.set_yticklabels(['2', '5', '10', '20', '50', '100'])
axe.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
GR = dict(fontsize=8, ha='left', va='center', arrowprops=dict(arrowstyle='->', lw=0.7, shrinkA=1, shrinkB=1))
GR['annotation_clip'] = False
tr = ('data', 'axes fraction')
axe.annotate('RANS, 0.9', xy=(0.9, 1.0), xycoords=tr, xytext=(0.2, 1.1), textcoords=tr, color='k', **dict(GR, ha='center'))
axe.annotate('K&B high FST, 1', xy=(1.0, 0.60), xycoords=tr, xytext=(1.2, 0.60), textcoords=tr, color='0.35', **GR)
axe.annotate('K&B low FST, $\\approx2$', xy=(2.0, 0.48), xycoords=tr, xytext=(2.2, 0.48), textcoords=tr, color='0.35', **GR)
axe.text(0.72, 0.90, 'solid: span-averaged\ndashed: centre plane', fontsize=8, color='k', ha='left', va='center', transform=axe.transAxes)
axe.text(0.72, 0.74, 'near film', fontsize=8, color='#0072B2', ha='left', va='center', transform=axe.transAxes)
axe.text(0.86, 0.74, 'far film', fontsize=8, color='#009E73', ha='left', va='center', transform=axe.transAxes)
axe.text(0.72, 0.66, r'own $u_\tau$ per region', fontsize=8, color='k', ha='left', va='center', transform=axe.transAxes)
axe.annotate('counter-gradient heat flux\n(faint, open markers)', xy=(1.6, 0.40), xycoords=tr, xytext=(2.4, 0.32),
             textcoords=tr, fontsize=8, ha='center', va='center', color='0.35',
             arrowprops=dict(arrowstyle='->', lw=0.7, shrinkA=1, shrinkB=1, color='0.35'))
fig.savefig(OUT, dpi=300)
print('saved', OUT)

# ---- printed results
print(f'grid: {CX * 0.02:.2f} D x {CZ * 0.02:.2f} D; empty cells {int((Cg["V"] <= 0).sum())}; well-defined share {K_RAW.mean():.2f} raw, {K1.mean():.2f} after island removal')
print(f'island removal: pieces below {MIN_AREA} D^2 (cell 0.004 D^2); removed {len(removed)}')
for val, area, xx, zz in removed:
    print(f'  {"well-defined island" if val else "undefined hole"}: area {area:.3f} D^2 at x/D {xx:.2f}, z/D {zz:.2f}')
print('wall units (span- and x-averaged mean wall shear, b010):')
for rn, lo, hi, _, _ in REG:
    w = WU[rn]
    print(f'  {rn} x/D {lo}-{hi}: <tau_w> {w["tau"]:.2f} Pa, <rho_w> {w["rho"]:.4f} kg/m3, nu_w {w["nu"]:.3e} m2/s, '
          f'u_tau {w["ut"]:.3f} m/s, mean C_f {w["cf"]:.5f}, l_nu = nu/u_tau {1e6 * w["nu"] / w["ut"]:.1f} um = {LS[rn]:.5f} D')
    print('   z+ for z/D 0.02 0.05 0.1 0.2 0.5 1.0:', ' '.join(f'{z / LS[rn]:.1f}' for z in (0.02, 0.05, 0.1, 0.2, 0.5, 1.0)))
    bad, badc, nu, al, nuc, alc = FLAG[rn]
    ps, _, pc, zp = prof[rn]
    for lab, bb, zz in (('span', bad, zs), ('centre', badc, zrow)):
        zf = zz[bb & (zz > 0) & (zz <= 0.5)]
        print(f'   flagged ({lab}), z/D <= 0.5: ' + (' '.join(f'{z:.2f}' for z in zf) if len(zf) else 'none') +
              f'  -> z+ {zf.min() / LS[rn]:.1f}-{zf.max() / LS[rn]:.1f}' if len(zf) else f'   flagged ({lab}): none')
    print('   z/D 0.02-0.10: nu_R=-A/B, alpha_R=-C/E (span):', ' '.join(f'{z:.2f}: {nu[i]:.2e}/{al[i]:.2e}' for i, z in zip(range(1, 6), zs[1:6])))
    print('   Pr_t span / centre at z+ 10, 30, 100:', ' | '.join(
        f'{zt}: {np.interp(zt, zp[1:], ps[1:]):.2f} / {np.interp(zt, zrow[1:] / LS[rn], pc[1:]):.2f}' for zt in (10, 30, 100)))
print('sign/fc check, all bins z/D >= 0.3:')
for rn, lo, hi in (('jet', 0, 3), ('near film', 3, 10), ('far film', 10, 20)):
    sel = (xs >= lo - 1e-9) & (xs < hi - 1e-9)
    zz = zs >= 0.3 - 1e-9
    t = [Q[k][np.ix_(sel, zz)].sum() for k in 'ABCE']
    print(f'  {rn} x/D {lo}-{hi}: Pr_t {ratio(*t):.3f}')
