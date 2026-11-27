"""fig_prt_classic_span: span-averaged version of fig_prt_classic (comparison file, centre-plane figure unchanged).

Pr_t(x, z) = (A/B) / (C/E) with the y-summed least-squares sums of the 193M 3DNS (stat3 b010), file
/mnt/Elements/DNS/analysis/asme2027_figs/n8_lsq/n8_lsq_c193b_b010.npz (0.02 D bins in x and z; written by
4dns_filmcooled/analysis/eps_prt_c193b/n8_sums.py):
    A = sum_y V rho u''w'' dU/dz,  B = sum_y V (dU/dz)^2,  C = sum_y V rho w''T'' dT/dz,  E = sum_y V (dT/dz)^2.
The ratio (A/B)/(C/E) needs no sign change and no density.
  (a) map over x/D -2..20, z/D 0..3. A, B, C, E are accumulated onto a regular grid of 0.1 D in x and 0.04 D in z
      (5 x 2 original bins, no empty cell) and the ratio is formed afterwards; nothing is interpolated or smoothed.
      Owner exception to the no-masking rule: cells outside the well-defined region are drawn in uniform light grey
      (#d9d9d9) and carry no Pr_t colour. The grey is set by exactly one criterion, the gradient test below
      (rms gradients sqrt(B/V) <= 0.05 U_inf/D or sqrt(E/V) <= 0.05 dT0/D), after the island removal, plus any cell
      with a non-finite ratio (none occurs). The counter-gradient and resolution flags belong to the profile
      computation only and do not enter the map. The covered share of the plotted area is printed ("grey share").
      Dashed line (#56B4E9): well-defined region, where the rms gradients of the cell, sqrt(B/V) and sqrt(E/V),
      exceed 0.05 U_inf/D (velocity) and 0.05 dT0/D (temperature), dT0 = 100.98 K; the threshold of the centre-plane
      figure applied to the span-rms gradient (V = sum of the grid-point volumes in the cell).
  The profile panel (Pr_t against z+) is no longer drawn; its values (span and centre plane, jackknife SE of the
  10 flow-through records from n8_lsq_c193b_b010_jk.npz, override with argv[1] or BATCH_FILE) are printed to stdout.
  Pr_t is formed per delete-one sample, SE = sqrt((n-1)/n sum (theta_b - mean)^2).
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
XSTART = float(os.environ.get('MAP_X0', '4'))      # left limit of the map in x/D (paper: 4; MAP_X0=-2 for the full range)
ZTOP = float(os.environ.get('MAP_ZTOP', '2.5' if XSTART >= 4 else '3'))   # top of the map in z/D (nothing is defined above z/D 2.2 for x/D 4-20)
# map data (2026-11-27): sums on 193M mesh bins (z: one bin per mesh level; x: 0.02 D upstream of x/D 10.25, one bin per
# mesh line downstream), merged in x to about 0.1 D; jk_* = ten delete-one jackknife samples. The coarse 5 x 2 map built
# from SUMS below is kept only for the old-against-new comparison printed at the end.
MESH = '/mnt/Elements/DNS/analysis/asme2027_figs/n8_lsq/n8_lsq_c193b_b010_meshbins_jk.npz'
XW = 0.1                                           # target x width of a map cell (D)
SUMS = '/mnt/Elements/DNS/analysis/asme2027_figs/n8_lsq/n8_lsq_c193b_b010_jk.npz'   # full-record sums, x/D -2..20 (identical to n8_lsq_c193b_b010.npz there)
JK = SUMS                           # single data path: switch SUMS to change the input file
BATCH_FILE = sys.argv[1] if len(sys.argv) > 1 else os.environ.get('BATCH_FILE', JK)
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
VMAX = 2.0
D = mc.D
NBOOT = 1000
CX, CZ = 5, 2                       # original bins per map cell in x and z (0.1 D x 0.04 D)
ZMAX_B = 1.0                     # upper z/D limit of panel (b)

set_mpl_params(10.0 / 22.0, grid_on=False)
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
# old coarse map, kept for the comparison
P_OLD, K_OLD, XE_OLD, ZE_OLD = P1.copy(), K1.copy(), xe.copy(), ze.copy()
K_RAW_OLD, removed_old = K_RAW, removed

# ---- new map from the mesh-bin file
M_ = np.load(MESH)
XED, ZED = M_['x_edges'].copy(), M_['z_edges'].copy()
ZED[0] = 0.0
nxn = len(XED) - 1
# merge columns: close a group when its width is nearer to XW than the width with the next column added (the 128-point
# bin at x/D 10.13-10.15 is merged with its neighbours this way)
starts, i0 = [0], 0
for i in range(1, nxn):
    w_cur, w_nxt = XED[i] - XED[i0], XED[i + 1] - XED[i0]
    if abs(w_cur - XW) <= abs(w_nxt - XW):
        starts.append(i)
        i0 = i
starts = np.array(starts)
XEM = np.concatenate([XED[starts], [XED[-1]]])                 # merged x edges
if XEM[-1] - XEM[-2] < 0.5 * XW:                               # short last cell: join the previous one
    XEM = np.delete(XEM, -2)
    starts = starts[:-1]


def merge(M):
    return np.add.reduceat(M, starts, axis=-2)


Cn = {k: merge(M_[k]) for k in 'ABCEV'}
Cn['N'] = merge(M_['N'])
JKn = {k: merge(M_['jk_' + k]) for k in 'ABCE'}
xe, ze = 0.5 * (XEM[1:] + XEM[:-1]), 0.5 * (ZED[1:] + ZED[:-1])
AREA = np.outer(np.diff(XEM), np.diff(ZED))
P1 = ratio(Cn['A'], Cn['B'], Cn['C'], Cn['E'])
with np.errstate(divide='ignore', invalid='ignore'):
    K1 = (np.sqrt(Cn['B'] / Cn['V']) > 0.05 * mc.U_INF / D) & (np.sqrt(Cn['E'] / Cn['V']) > 0.05 * mc.DT0 / D)
K_RAW = K1.copy()
removed = []
for val in (True, False):                                       # same island removal as above, with true cell areas
    lab, nl = ndimage.label(K1 == val)
    for i in range(1, nl + 1):
        m_ = lab == i
        area = AREA[m_].sum()
        if area < MIN_AREA:
            ii, jj = np.nonzero(m_)
            removed.append((val, area, xe[ii].mean(), ze[jj].mean()))
            K1[m_] = not val
# jackknife standard error of Pr_t per cell (delete-one, n = 10) and the two resolution criteria
TH = ratio(JKn['A'], JKn['B'], JKn['C'], JKn['E'])             # (10, nx, nz)
NJ = TH.shape[0]
with np.errstate(invalid='ignore'):
    SE = np.sqrt((NJ - 1) / NJ * ((TH - TH.mean(0)) ** 2).sum(0))
    UNRES_SE = ~(np.abs(P1) >= 2 * SE)                         # also true where SE is not finite
SIGN = np.zeros(P1.shape, bool)
for k in 'BC':                                                  # denominators of (A/B)/(C/E); E > 0 by construction
    s_ = np.concatenate([JKn[k], Cn[k][None]])
    SIGN |= ~((s_ > 0).all(0) | (s_ < 0).all(0))
XC, ZC = np.meshgrid(xe, ze, indexing='ij')
XEg, ZEg = np.meshgrid(XEM, ZED, indexing='ij')

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


# +-2 SE band: delete-one jackknife over the 10 flow-through records (jk_A.. arrays [10, nx, nz])
band, JKSE, JKAL = {}, {}, {}
if BATCH_FILE and os.path.exists(BATCH_FILE):
    bz = np.load(BATCH_FILE)
    for rn, lo, hi, _, _ in REG:
        sel = (xs >= lo - 1e-9) & (xs < hi - 1e-9)
        sb = [bz['jk_' + k][:, sel].sum(1) for k in 'ABCE']             # 4 x (nb, nz), delete-one sums
        nb = sb[0].shape[0]
        th = ratio(*sb)                                                 # (nb, nz) Pr_t per delete-one sample
        se = np.sqrt((nb - 1) / nb * ((th - th.mean(0)) ** 2).sum(0))
        al_b = -sb[2] / sb[3]                                           # alpha_R per delete-one sample
        al_se = np.sqrt((nb - 1) / nb * ((al_b - al_b.mean(0)) ** 2).sum(0))
        nu_b = -sb[0] / sb[1]                                           # nu_R per delete-one sample
        nu_se = np.sqrt((nb - 1) / nb * ((nu_b - nu_b.mean(0)) ** 2).sum(0))
        band[rn] = se
        JKSE[rn] = th
        JKAL[rn] = (al_b, al_se, nu_b, nu_se)

# Flow-field map at full text width (owner 2026-11-27): 7.5 in wide (figure*, width=\textwidth, scale 1), 10 pt text.
# MAP_ZX = z stretch factor of the axes (1 = true aspect: equal scales in x and z; at most 2); MAP_OUT overrides the output file.
ZX = float(os.environ.get('MAP_ZX', '2'))      # default 2: the stretched version reads better (annotations clear of the boundary)
OUT = os.environ.get('MAP_OUT', OUT)
W = 7.5
AX_L, AX_W, CB_GAP, CB_W = 0.55, 6.3, 0.12, 0.11            # inches: left margin, axes width, gap to the colour bar, colour bar width
AX_H = AX_W * ZTOP / (20.0 - XSTART) * ZX                      # x/D XSTART..20 by z/D 0..ZTOP at equal scales, times the stretch
AX_B, AX_T = 0.43, 0.08
H = AX_B + AX_H + AX_T
fig = plt.figure(figsize=(W, H))
kw = dict(cmap='inferno', vmin=0.0, vmax=VMAX, shading='auto', rasterized=True)
ckw = dict(levels=[0.5], colors='#56B4E9', linewidths=0.9, linestyles='--')
ax1 = fig.add_axes([AX_L / W, AX_B / H, AX_W / W, AX_H / H])
cax = fig.add_axes([(AX_L + AX_W + CB_GAP) / W, AX_B / H, CB_W / W, AX_H / H])
print(f'figsize {W} x {H:.2f} in, axes {AX_W:.2f} x {AX_H:.2f} in, z stretch {ZX:g}')

GREY = '#d9d9d9'
UNDEF = ~K1 | ~np.isfinite(P1)                                  # present criterion: gradient test (after island removal), non-finite
GREY_MASK = UNDEF | UNRES_SE | SIGN                             # plus: |Pr_t| < 2 SE, or a denominator sum changes sign among the samples
GREY_FRAC = (AREA * GREY_MASK).sum() / AREA.sum()               # cells differ in size: area-weighted
m = ax1.pcolormesh(XEg, ZEg, np.ma.masked_array(P1, mask=GREY_MASK), **kw)
ax1.contour(XC, ZC, K1.astype(float), **ckw)
ax1.set_facecolor(GREY)
ax1.annotate(r'$\mathrm{Pr}_t$ not defined', xy=(10.5, 2.2 if ZTOP < 3 else 2.35), xytext=(19.6, 2.3 if ZTOP < 3 else 2.6),
             ha='right', va='center', color='black',
             arrowprops=dict(arrowstyle='->', color='black', lw=0.8, shrinkA=1, shrinkB=1))
ax1.set_ylabel(r'$z/D$', labelpad=1)
ax1.set_yticks(list(range(int(ZTOP) + 1)))
ax1.set_ylim(0, ZTOP)
ax1.set_xlim(XSTART, 20)
ax1.set_xticks([t for t in (0, 5, 10, 15, 20) if t >= XSTART])
ax1.set_xlabel(r'$x/D$', labelpad=1)
ax1.annotate('boundary of the\nwell-defined region', xy=(15.0, 1.95), xytext=(15.0, 1.15),
             ha='center', va='center', color='white',
             arrowprops=dict(arrowstyle='->', color='white', lw=0.8, shrinkA=1, shrinkB=1))
cb = fig.colorbar(m, cax=cax, orientation='vertical', extend='both')
cb.set_ticks([0, 0.5, 0.9, 1.5, 2.0])
cb.set_ticklabels(['0', '0.5', '0.9', '1.5', '2'])
cb.ax.axhline(0.9, color='white', lw=1.4)
cb.set_label(r'$\mathrm{Pr}_t$', labelpad=2)

# ---- profiles (not plotted): Pr_t against z+, each region with its own u_tau; values are printed below
prof = {}
zp0 = zs > 0                                   # the wall bin z = 0 has no z+
FLAG = {}
ZLO = {}
for rn, lo, hi, col, ls in REG:
    ps, bad, nu, al = span_profile(lo, hi)
    pc, badc, nuc, alc = centre_profile(lo, hi)
    zp = zs / LS[rn]
    prof[rn] = (ps, bad, pc, zp)
    FLAG[rn] = (bad, badc, nu, al, nuc, alc)
    if rn in band:
        al_b, al_se, nu_b, nu_se = JKAL[rn]
        res = ~bad & (np.abs(al) > 2 * al_se) & (np.abs(nu) > 2 * nu_se)      # resolved bins
        res = res | (zs > ZMAX_B + 1e-9)
        up = np.cumprod(res[::-1])[::-1].astype(bool)                 # resolved from here up to the top
        i0 = int(np.nonzero(up & zp0)[0].min())
        ZLO[rn] = (zp[i0], ps[i0])
if ZLO:
    print('lowest resolved z+:', {k: round(v[0], 1) for k, v in ZLO.items()})
fig.savefig(OUT, dpi=600)
print('saved', OUT)
print(f'grey share (not defined or not resolved, x/D -2..20, z/D 0..3): {GREY_FRAC:.3f} of the plotted area')

# ---- printed results
print(f'grid: {len(xe)} x cells ({np.diff(XEM).min():.3f}-{np.diff(XEM).max():.3f} D, median {np.median(np.diff(XEM)):.3f}) by {len(ze)} z levels; '
      f'empty cells {int((Cn["V"] <= 0).sum())}; min points per cell {int(Cn["N"].min())}; well-defined area share '
      f'{(AREA * K_RAW).sum() / AREA.sum():.3f} raw, {(AREA * K1).sum() / AREA.sum():.3f} after island removal')

# ---- grey shares by criterion (area-weighted, x/D -2..20, z/D 0..3) and in the far-film near-wall strip
def shares(sel, label):
    a = (AREA * sel).sum()
    f = lambda m: (AREA * m * sel).sum() / a
    print(f'  {label}: area {a:.3f} D^2; undefined (gradient test) {f(UNDEF):.3f}; |Pr_t| < 2 SE {f(UNRES_SE):.3f}; '
          f'denominator sign change {f(SIGN):.3f}; union {f(GREY_MASK):.3f}; extra to undefined: SE-only {f(~UNDEF & UNRES_SE & ~SIGN):.3f}, '
          f'sign-only {f(~UNDEF & ~UNRES_SE & SIGN):.3f}, either {f(~UNDEF & (UNRES_SE | SIGN)):.3f}')


XX, ZZ = np.meshgrid(xe, ze, indexing='ij')
print('grey shares of the drawn area:')
shares(np.ones_like(P1, bool), 'all, x/D -2..20')
shares((XX >= XSTART), f'drawn range x/D {XSTART:g}..20')
shares((XX >= 10) & (XX < 20) & (ZZ < 0.2), 'far-film near-wall strip x/D 10-20, z/D < 0.2')
shares((XX >= 4) & (XX < 10) & (ZZ < 0.2), 'near-film near-wall strip x/D 4-10, z/D < 0.2')
shares((XX >= 10) & (XX < 20) & (ZZ >= 0.2), 'far film x/D 10-20, z/D >= 0.2')
shares((XX >= 4) & (XX < 10) & (ZZ >= 0.2), 'near film x/D 4-10, z/D >= 0.2')
print('median SE of Pr_t, strip x/D 10-20, z/D < 0.2, defined cells:',
      f'{np.nanmedian(SE[(XX >= 10) & (XX < 20) & (ZZ < 0.2) & ~UNDEF]):.3f}')


# ---- old (5 x 2 merged, 0.1 D x 0.04 D) against new
def wq(v, w, q):
    o = np.argsort(v)
    v, w = v[o], w[o]
    c = (np.cumsum(w) - 0.5 * w) / w.sum()
    return np.interp(q, c, v)


XO, ZO = np.meshgrid(XE_OLD, ZE_OLD, indexing='ij')
print('Pr_t median [10-90 %] over defined cells (old: K_OLD; new: area-weighted, defined = gradient test; resolved = also SE and sign):')
for lab, xl, xh, zl, zh in (('x/D 4-10, z/D >= 0.3', 4, 10, 0.3, 9), ('x/D 10-20, z/D >= 0.3', 10, 20, 0.3, 9),
                            ('strip x/D 4-20, z/D < 0.2', 4, 20, 0.0, 0.2)):
    so = (XO >= xl) & (XO < xh) & (ZO >= zl) & (ZO < zh) & K_OLD & np.isfinite(P_OLD)
    sn = (XX >= xl) & (XX < xh) & (ZZ >= zl) & (ZZ < zh) & ~UNDEF
    sr = sn & ~UNRES_SE & ~SIGN
    out = []
    for nm, s_, P_, w_ in (('old', so, P_OLD, np.ones_like(P_OLD)), ('new defined', sn, P1, AREA), ('new resolved', sr, P1, AREA)):
        q = wq(P_[s_], w_[s_], [0.1, 0.5, 0.9])
        out.append(f'{nm} {q[1]:.3f} [{q[0]:.3f}, {q[2]:.3f}] (n {int(s_.sum())})')
    print(f'  {lab}: ' + '; '.join(out))
print('volume-summed Pr_t, z/D >= 0.3, all native mesh bins (compare 0.300, 0.225, 0.454 of the earlier map):')
for rn, lo, hi in (('jet', 0, 3), ('near film', 3, 10), ('far film', 10, 20)):
    sx = (M_['x_D'] >= lo) & (M_['x_D'] < hi)
    sz = M_['z_D'] >= 0.3
    t = [M_[k][np.ix_(sx, sz)].sum() for k in 'ABCE']
    print(f'  {rn} x/D {lo}-{hi}: Pr_t {ratio(*t):.3f}')
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
if band:
    print('jackknife (delete-one, n = 10) Pr_t span-averaged, SE; band drawn is +-2 SE')
    for rn, lo, hi, _, _ in REG:
        ps, bad, nu, al = span_profile(lo, hi)
        se = band[rn]
        zp = zs / LS[rn]
        print(f'  {rn} x/D {lo}-{hi}: z/D, z+, Pr_t, SE')
        for zd in (0.1, 0.2, 0.3, 0.5):
            i = int(round(zd / 0.02))
            print(f'    z/D {zs[i]:.2f}  z+ {zp[i]:6.1f}  Pr_t {ps[i]:.3f}  SE {se[i]:.3f}  flagged {bool(bad[i])}')
        for zt in (30, 100):
            i = int(np.argmin(np.abs(zp[1:] - zt))) + 1
            print(f'    z/D {zs[i]:.2f}  z+ {zp[i]:6.1f}  Pr_t {ps[i]:.3f}  SE {se[i]:.3f}  flagged {bool(bad[i])}')
    rn = 'far film'
    al_b, al_se = JKAL[rn][:2]
    _, bad, nu, al = span_profile(10, 20)
    zp = zs / LS[rn]
    print('counter-gradient, far film (alpha_R = -C/E, full record; negative in n of 10 delete-one samples):')
    for i in np.nonzero((zp >= 2.5) & (zp <= 12) & (zs > 0))[0]:
        print(f'    z/D {zs[i]:.2f}  z+ {zp[i]:5.2f}  alpha_R {al[i]:.3e}  SE {al_se[i]:.2e}  n(alpha_R<0) {int((al_b[:, i] < 0).sum())}/10  flagged {bool(bad[i])}')
print('sign/fc check, all bins z/D >= 0.3:')
for rn, lo, hi in (('jet', 0, 3), ('near film', 3, 10), ('far film', 10, 20)):
    sel = (xs >= lo - 1e-9) & (xs < hi - 1e-9)
    zz = zs >= 0.3 - 1e-9
    t = [Q[k][np.ix_(sel, zz)].sum() for k in 'ABCE']
    print(f'  {rn} x/D {lo}-{hi}: Pr_t {ratio(*t):.3f}')
