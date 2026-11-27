"""F14: span-mean maps of the factors of the eddy viscosity (RANS 50M SST against the 193M reference b010) and an exact factor split.

Usage: python3 fig_meff_maps.py
Output: /home/sk999/Documents/ASME_2027/figures/fig_meff_maps.png, data/meff_maps_table.txt, data/near_wall_boundary*.txt,
        /home/sk999/Documents/ASME_2027/planning/review/24_factor_maps_raw.png (review sheet, raw factor fields)

mu = rho k^(1/2) L cos(theta) at every point. Lengths: 193M L_eff = k^(1/2)|b|/|S*| (L of the point file, mu_eff = rho k^(1/2) L_eff cos);
RANS L_R = mu_t/(rho k^(1/2)) (cos = 1 by construction).
Maps (x/D -2..20, z/D 0..2), ratios RANS/193M of plain span means over the full pitch, one log10 bwr scale (0.1 to 10, white = 1), so that
(a) k^(1/2), (b) length, (c) 1/cos(theta) (cos itself is signed and would need a different scale; 1/cos on the common log scale reads
better because the factor ratios then multiply to (d)) and the density ratio (about 1, not drawn) multiply to (d) the ratio of the span-mean mu.
Thin black contour = ratio 1. (a)-(c): span means of the factors computed from the point file (the meff_span files carry only mu,
S_mean, S_frac_gt005); complete span only for x/D >= 2.8, blank (jet region) elsewhere. (d): span means of mu from the meff_span files.
Grey (owner exception to the no-masking rule, as for the Pr_t map): above the strain boundary z_s(x), where the span-mean 193M
|S*| D/U_inf < 0.05 (mu_eff = -a:S*/(2|S*|^2) undefined without mean strain), and, in (a)-(d), where the ratio is not defined (mean <= 0).
(e) bars, WEIGHTING = 'mask' (default): strain acts only as a mask (point-wise 193M |S*| D/U_inf > 0.05 kept), then plain
volume-weighted statistics: total = ln(sum V mu_t / sum V mu_eff) and its exact Shapley split into rho, k^(1/2), L_eff, cos (four
players, no strain-weight player), x/D 3-10 near film and 10-20 far film, z <= 2 D, wall nodes excluded.
WEIGHTING = 'strain2': previous strain-squared weighted regional least squares with five players (and BOUNDARY layers).

BOUNDARY (switch below) selects the near-wall / far-wall boundary used for the table layers and the lines:
  'zx'        z_x(x), where the two span means are equal.
  'strain005' z_s(x), lowest height above the wall where the 193M span-mean |S*| D/U_inf falls below 0.05 (default).
Source of the span means of mu: meff_span files (RANS50 + C193b) when present, otherwise span_fields.npz for both.
"""
BOUNDARY = 'strain005'      # 'strain005' (default) or 'zx'
WEIGHTING = 'mask'          # 'mask' (default) or 'strain2' (previous)
MASK = '3d'                 # '3d' (default): point-wise strain mask applied before the span mean (conditional span means, point file);
                            # '2d' (previous): span means of unconditional fields, grey above z_s of the span-mean strain
import os
ALT = bool(os.environ.get('F14_ALT'))          # alternative version for the owner: f_min = 0.2, x axis from x/D 4 (figure only, no tables)
SRC = os.environ.get('F14_SRC', 'pointfile')     # 'pointfile' (default, method of picture 25); 'span' (alternative): full-pitch conditional means of the meff_span files (c_*_mean), RANS interpolated onto
                                            # the 193M grid; 'pointfile': previous source (point file pts.npz, reduced set of points)
FMIN = float(os.environ.get('F14_FMIN', 0.2))   # 3d: a map cell is drawn only where the valid fraction f of the span points is >= FMIN
XMIN = float(os.environ.get('F14_XMIN', 4.0))   # x/D < 4 is left out of the maps, z_x and the split (owner 2026-11-27; was 3)
# MANUAL MASK OF THE JET EDGE (drawing only; z_x, the split and the table do not use it): the cells above the straight line from
# (JET_MASK_X0, JET_MASK_Z0) to (JET_MASK_X1, JET_MASK_Z1) are drawn grey in all four panels. The line joins the existing grey area at
# x/D JET_MASK_X1 and hides the saturated patch (ratio >= 3 in panels c and d) at x/D 3-4.05, z/D 1.23-1.53.
JET_MASK_X0, JET_MASK_Z0 = 3.0, 1.15
JET_MASK_X1, JET_MASK_Z1 = 4.25, 1.41
LAYOUT = 'row'            # 'row': maps (a)-(c) in one row, (d) full width, (e) beneath; 'stack': four maps stacked
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Patch

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from functions import set_mpl_params
import _meff_factor_cond as FC
import _meff_factor_span as SF
import _meff_pointwise_data as M
import _meff_span_split as SS2
import _meff_split_nearfar as SP
import _meff_strain_span as SS

OUT = (os.environ.get('F14_OUT') or (f'/home/sk999/Documents/ASME_2027/planning/review/26_meff_maps_fullpitch_fmin{str(FMIN).replace(".", "").rstrip("0") if FMIN != 0.2 else "02"}.png' if ALT
       else '/home/sk999/Documents/ASME_2027/figures/fig_meff_maps.png'))
RAW = '/home/sk999/Documents/ASME_2027/planning/review/24_factor_maps_raw.png'
TAB = '/home/sk999/Documents/ASME_2027/figures/code/data/meff_maps_table.txt'
SPAN_DIR = os.environ.get('F14_SPAN_DIR', '/mnt/Elements/DNS/analysis/asme2027_figs/meff_span/')   # env override only for tests
FALLBACK = M.SCR + 'span_fields.npz'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
XT0, ZT0 = np.array([3, 5, 10, 15, 20.]), np.array([0.62, 0.61, 0.72, 0.87, 0.94])
ZMAX = 2.0
UINF, DD = M.UINF, M.D


def load_source():
    f_r, f_d = SPAN_DIR + 'span_RANS50.npz', SPAN_DIR + 'span_C193b.npz'
    if not (os.path.exists(f_r) and os.path.exists(f_d) and 'S_mean' in np.load(f_d).files):
        v1r, v1d = SPAN_DIR + 'v1_no_strain/span_RANS50.npz', SPAN_DIR + 'v1_no_strain/span_C193b.npz'
        if os.path.exists(v1r) and os.path.exists(v1d):
            print('WARNING: span_C193b.npz with S_mean not present; using the previous files in v1_no_strain/ (no strain field)')
            f_r, f_d = v1r, v1d
    if os.path.exists(f_r) and os.path.exists(f_d):
        r, d = np.load(f_r), np.load(f_d)
        from scipy.interpolate import RegularGridInterpolator as RGI
        mu = float(d['mu_ref'])
        # the two files have different bins (RANS: 0.04 D in x, 110 z levels; 193M: 0.02 D, 163 z levels): both are interpolated
        # linearly (bin centres) onto one regular grid, x/D step 0.04, z/D step 0.01
        xc = np.arange(-2 + 0.02, 20 + 1e-9, 0.04)
        zc = np.arange(0, 2.0 + 1e-9, 0.01)
        ze = np.r_[0.0, 0.5 * (zc[1:] + zc[:-1]), zc[-1] + 0.005]
        xe = np.r_[xc - 0.02, xc[-1] + 0.02]

        def regrid(f):
            ip = RGI((f['z_D'], f['x_D']), (f['mean'].T / mu), bounds_error=False, fill_value=None, method='linear')
            q = np.stack(np.meshgrid(np.clip(zc, f['z_D'][0], f['z_D'][-1]), np.clip(xc, f['x_D'][0], f['x_D'][-1]), indexing='ij'), -1)
            return ip(q)

        Rm, Dm = regrid(r), regrid(d)                                    # [nz, nx]
        good = np.isfinite(Rm) & np.isfinite(Dm)
        Smean = Sfrac = None
        if 'S_mean' in d.files:                                          # strain of the same file as the 193M maps
            def regrid_k(f, key):
                ip = RGI((f['z_D'], f['x_D']), f[key].T, bounds_error=False, fill_value=None, method='linear')
                return ip(np.stack(np.meshgrid(np.clip(zc, f['z_D'][0], f['z_D'][-1]), np.clip(xc, f['x_D'][0], f['x_D'][-1]),
                                               indexing='ij'), -1))
            Smean, Sfrac = regrid_k(d, 'S_mean'), regrid_k(d, 'S_frac_gt005')
        src = ('meff_span: span_RANS50.npz (RANS 50M mu_t) and span_C193b.npz (193M b010 mu_eff); plain mean over all grid points of the '
               'full pitch, wall nodes included (mu = 0 at the wall row); the two files (different bins) interpolated linearly onto one '
               'grid, x/D step 0.04, z/D step 0.01; z/D 0-3 available (maps drawn to z/D 2)')
        short = 'meff_span (RANS50 + C193b)' + (' with S_mean' if Smean is not None else ', no strain field (v1_no_strain)' if 'v1_no_strain' in f_d else '')
        return dict(Rm=Rm, Dm=Dm, xe=xe, ze=ze, xc=xc, zc=zc, good=good, src=src, short=short, incomplete=None, Smean=Smean, Sfrac=Sfrac)
    f = np.load(FALLBACK)
    R, Dd, x, z = f['R'], f['D'], f['x'], f['z']
    nR, nD = np.isfinite(R).sum(2), np.isfinite(Dd).sum(2)
    good = (nR >= 119) & (nD >= 119)                                    # complete span (>= 119 of 121 points) in both
    with np.errstate(all='ignore'):
        Rm, Dm = np.nanmean(R, axis=2), np.nanmean(Dd, axis=2)
    Rm, Dm = np.where(good, Rm, np.nan), np.where(good, Dm, np.nan)
    xe = np.r_[x - 0.025, x[-1] + 0.025]
    ze = np.r_[z - 0.005, z[-1] + 0.005]
    ze[0] = 0.0
    cx = x[good[60]]
    src = ('FALLBACK span_fields.npz (both rungs): point-file span means interpolated to a regular grid, plain mean over the points '
           f'of the span line, complete span only (>= 119 of 121 points); x/D < {cx.min():.2f} blank; z/D <= 2 only; '
           'meff_span/span_C193b.npz not present')
    return dict(Rm=Rm, Dm=Dm, xe=xe, ze=ze, xc=x, zc=z, good=good, src=src, short='fallback span_fields.npz (both rungs)',
                incomplete=cx.min(), Smean=None, Sfrac=None)


S = load_source()
print('SOURCE:', S['short'])
print('  ', S['src'])
Rm, Dm, xc, zc, xe, ze, good = S['Rm'], S['Dm'], S['xc'], S['zc'], S['xe'], S['ze'], S['good']
nz, nx = Rm.shape
dx, dz = np.diff(xe), np.diff(ze)
XC, ZC = np.meshgrid(xc, zc)
AREA = np.outer(dz, dx)

# low-strain mask (owner exception to the no-masking rule, as for the Pr_t map): cells where the span-mean strain magnitude of the
# 193M, |S*| D/U_inf, is below S0 are drawn grey (mu_eff = -a:S*/(2|S*|^2) is not defined without mean strain).
S0 = 0.05
Smean, Sfrac = S['Smean'], S['Sfrac']
if Smean is not None:
    print('STRAIN SOURCE: S_mean / S_frac_gt005 of span_C193b.npz (same file as the maps)')
else:
    xs_, zs_, Smean, Sfrac = SS.get(S0)
    print('STRAIN SOURCE: point file (span_fields grid)')
    if S['incomplete'] is None:
        print('WARNING: the meff_span 193M file carries no |S*| field; the low-strain mask comes from the point file '
              '(span_fields grid) and is mapped to the meff_span bins by nearest neighbour. Request a span-mean |S*| field.')
        ix = np.abs(xc[:, None] - xs_[None, :]).argmin(1)
        iz = np.abs(zc[:, None] - zs_[None, :]).argmin(1)
        Smean, Sfrac = Smean[np.ix_(iz, ix)], Sfrac[np.ix_(iz, ix)]
LOW = np.isfinite(Smean) & (Smean < S0)          # criterion: span-mean |S*| D/U_inf < 0.05
LOWALT = np.isfinite(Sfrac) & (Sfrac < 0.5)      # alternative: fraction of span points with |S*| D/U_inf > 0.05 below one half
Rd, Dd_ = np.where(LOW, np.nan, Rm), np.where(LOW, np.nan, Dm)    # defined cells only
print(f'low-strain mask: span-mean |S*|D/U_inf < {S0}: {100 * LOW[good & (ZC > 0) & (ZC <= ZMAX) & (XC >= -2) & (XC <= 20)].mean():.1f} % of cells with data')


def crossing(col_r, col_d, z):
    """lowest upward crossing of RANS - 193M (below: RANS lower) with linear interpolation; z > 0 rows only."""
    df = col_r - col_d
    for i in range(1, len(z) - 1):
        a, b = df[i], df[i + 1]
        if np.isfinite(a) and np.isfinite(b) and a < 0 <= b:
            return z[i] + (0 - a) / (b - a) * (z[i + 1] - z[i])
    return np.nan


zrows = zc <= ZMAX + 1e-9


def zx_line(R_, D_):
    zc_ = np.array([crossing(R_[zrows, j], D_[zrows, j], zc[zrows]) for j in range(nx)])
    zs = np.full(nx, np.nan)          # smoothed z_x(x): running median over +-0.25 in x/D (the raw crossing is ragged in x)
    for j in range(nx):
        w = (np.abs(xc - xc[j]) <= 0.25) & np.isfinite(zc_)
        if w.sum() >= 3:
            zs[j] = np.median(zc_[w])
    return zs


zxs_old = zx_line(Rm, Dm)
okx_old = np.isfinite(zxs_old)
ZXT_OLD = np.interp(XT0, xc[okx_old], zxs_old[okx_old])
zxs = zx_line(Rd, Dd_)                    # defined cells only
okx = np.isfinite(zxs)
ZXT = np.interp(XT0, xc[okx], zxs[okx])
print('z_x(x) at x/D 3, 5, 10, 15, 20, all cells (old):   ' + ', '.join(f'{v:.3f}' for v in ZXT_OLD))
print('z_x(x) at x/D 3, 5, 10, 15, 20, defined cells (new): ' + ', '.join(f'{v:.3f}' for v in ZXT))
print('  split table z_x:                  ' + ', '.join(f'{v:.2f}' for v in ZT0))
dmax = np.abs(ZXT - ZT0).max()
RERUN = dmax > 0.03
print(f'  max |difference| = {dmax:.3f} D; split rerun with the new z_x: {RERUN}')
ZSPLIT = ZXT if RERUN else ZT0
zx_x = lambda xq: np.interp(xq, xc[okx], zxs[okx])  # noqa: E731
ZXC = zx_x(xc)

# z_s(x): lowest height above the wall where the span-mean 193M |S*| D/U_inf falls below S0 (linear interpolation), x/D >= 3,
# running median over +-0.25 in x/D
zs_raw = np.full(nx, np.nan)
for j in range(nx):
    col = Smean[:, j]
    for i in range(1, len(zc)):
        if zc[i] > ZMAX + 1e-9:
            break
        if np.isfinite(col[i]) and np.isfinite(col[i - 1]) and col[i - 1] >= S0 > col[i]:
            zs_raw[j] = zc[i - 1] + (col[i - 1] - S0) / (col[i - 1] - col[i]) * (zc[i] - zc[i - 1])
            break
zss = np.full(nx, np.nan)
for j in range(nx):
    w = (np.abs(xc - xc[j]) <= 0.25) & np.isfinite(zs_raw) & (xc >= 2.75)
    if w.sum() >= 3 and xc[j] >= 3:
        zss[j] = np.median(zs_raw[w])
oks = np.isfinite(zss)
ZST = np.interp(XT0, xc[oks], zss[oks])
ZSC = np.interp(xc, xc[oks], zss[oks])
print('z_s(x) at x/D 3, 5, 10, 15, 20 (|S*|D/U = 0.05): ' + ', '.join(f'{v:.3f}' for v in ZST))
print('z_x(x) at x/D 3, 5, 10, 15, 20 (RANS = 193M):    ' + ', '.join(f'{v:.3f}' for v in ZXT))
STRAIN = BOUNDARY == 'strain005'
MASKW = WEIGHTING == 'mask'
print(f'BOUNDARY = {BOUNDARY}, WEIGHTING = {WEIGHTING}')
INF = np.full(nx, np.inf)
ZERO = np.full(nx, -np.inf)
if STRAIN:
    LAY = {'near wall': (ZERO, ZSC), 'near wall lower': (ZERO, ZXC), 'near wall upper': (ZXC, ZSC), 'far wall': (ZSC, INF)}
    xs3 = xc[(xc >= 3) & oks]
    TB = lambda z: (xs3, np.interp(xs3, xc, z))  # noqa: E731
    SPLIT_LAYERS = {'near wall': (None, TB(ZSC)), 'near wall lower': (None, TB(ZXC)), 'near wall upper': (TB(ZXC), TB(ZSC)),
                    'far wall': (TB(ZSC), None)}
    KEYS = [('near film', 'near wall'), ('far film', 'near wall')]
    TKEYS = [(r_, l_) for r_ in ('near film', 'far film') for l_ in ('near wall', 'near wall lower', 'near wall upper', 'far wall')]
else:
    LAY = {'near wall': (ZERO, ZXC), 'far wall': (ZXC, INF)}
    SPLIT_LAYERS = None
    KEYS = [('near film', 'near wall'), ('near film', 'far wall'), ('far film', 'near wall'), ('far film', 'far wall')]
    TKEYS = KEYS

# band ratios (ratio of band means, bin-area weights; cells of the maps, wall row excluded, complete data)
BANDS = {}
GREY = {}
ZXC_OLD = np.interp(xc, xc[okx_old], zxs_old[okx_old])
for rk, lk in TKEYS:
    lo, hi = (3, 10) if rk == 'near film' else (10, 20.0001)

    def layer(zxc):
        if lk == 'near wall' and not STRAIN:
            lay = ZC < zxc[None, :]
        elif lk == 'far wall' and not STRAIN:
            lay = ZC >= zxc[None, :]
        else:
            lay = (ZC >= LAY[lk][0][None, :]) & (ZC < LAY[lk][1][None, :])
        return good & (XC >= lo) & (XC < hi) & lay & (ZC > 0) & (ZC <= ZMAX) & np.isfinite(Rm) & np.isfinite(Dm)
    m0 = layer(ZXC_OLD)
    a0, b0 = (Rm[m0] * AREA[m0]).sum(), (Dm[m0] * AREA[m0]).sum()
    m = layer(ZXC) & ~LOW                                           # defined cells only
    a, b = (Rm[m] * AREA[m]).sum(), (Dm[m] * AREA[m]).sum()
    mm_ = layer(ZXC)
    gA = AREA[mm_]
    nonpos = (Rm[mm_] <= 0) | (Dm[mm_] <= 0)
    GREY[(rk, lk)] = dict(low=(gA * LOW[mm_]).sum() / gA.sum(), merged=(gA * (LOW[mm_] | nonpos)).sum() / gA.sum(),
                          alt=(gA * LOWALT[mm_]).sum() / gA.sum(), both=(gA * (LOW[mm_] & LOWALT[mm_])).sum() / gA.sum())
    BANDS[(rk, lk)] = dict(ratio=a / b, mR=a / AREA[m].sum(), mD=b / AREA[m].sum(), n=int(m.sum()), old=a0 / b0)
    print(f'band {rk:9s} {lk:16s}: ratio all cells (old) {a0 / b0:.2f}; defined cells only (new) {a / b:.2f}')
areaw = (good & (ZC > 0) & (ZC <= ZMAX) & (XC >= XMIN) & (XC <= 20))

# strain-squared weighted split (previous definition, kept for the record)
from time import time
t0 = time()
SPL = SP.compute(XT0, ZSPLIT, layers=SPLIT_LAYERS)
print(f'strain-weighted split computed in {time() - t0:.0f} s')
# strain as a mask only (new primary definition)
MK = SP.compute_mask(S0, ZMAX)
for k, r in MK.items():
    print(f'mask split {k[0]:9s} {k[1]:7s}: total {r["tot"]:+.3f} rho/k/L/cos ' + ' '.join(f'{v:+.3f}' for v in r['p4']) +
          f' | residual {r["res"]:+.1e} | kept {r["kept"]:.3f} top {r["top"]:.2f} fneg {r["fneg"]:.3f}')

# ---------------------------------------------------------------- factor span means (point file)
xf, zf, FF, nf = SF.get()
pf = np.load(FALLBACK)
with np.errstate(all='ignore'):
    Rpf, Dpf = np.nanmean(pf['R'], axis=2), np.nanmean(pf['D'], axis=2)
avail = nf >= 119                                   # complete span
XF, ZF = np.meshgrid(xf, zf)
with np.errstate(all='ignore'):
    RK = FF['k12R'] / FF['k12D']
    RLn = FF['LR'] / FF['LD']
    RC = 1.0 / FF['cD']
    RRHO = FF['rhoR'] / FF['rhoD']
    RMU = Rpf / Dpf
    PROD = RRHO * RK * RLn * RC
FSRC = 'point file (pts.npz): span means interpolated as span_compute.py, complete span only (x/D >= 2.8), cached span_factors.npz'
_fr, _fd = SPAN_DIR + 'span_RANS50.npz', SPAN_DIR + 'span_C193b.npz'
if (MASK == '2d' and S['incomplete'] is None and os.path.exists(_fr) and os.path.exists(_fd) and 'sqrtk_mean' in np.load(_fr).files
        and 'sqrtk_mean' in np.load(_fd).files):
    from scipy.interpolate import RegularGridInterpolator as RGI

    def _rg(f, key):
        ip = RGI((f['z_D'], f['x_D']), f[key].T, bounds_error=False, fill_value=None, method='linear')
        return ip(np.stack(np.meshgrid(np.clip(zc, f['z_D'][0], f['z_D'][-1]), np.clip(xc, f['x_D'][0], f['x_D'][-1]), indexing='ij'), -1))
    _r, _d = np.load(_fr), np.load(_fd)
    with np.errstate(all='ignore'):
        RK = _rg(_r, 'sqrtk_mean') / _rg(_d, 'sqrtk_mean')
        RLn = _rg(_r, 'Leff_mean') / _rg(_d, 'Leff_mean')
        RC = 1.0 / _rg(_d, 'cos_mean')                       # RANS cos = 1 by construction
        RRHO = _rg(_r, 'rho_mean') / _rg(_d, 'rho_mean')
        RMU = Rm / Dm
        PROD = RRHO * RK * RLn * RC
    xf, zf, XF, ZF = xc, zc, XC, ZC
    avail = np.isfinite(Rm) & np.isfinite(Dm)
    FF = dict(cD=1.0 / RC, k12R=np.nan * RK, k12D=np.nan * RK, LR=np.nan * RK, LD=np.nan * RK, rhoR=RRHO, rhoD=np.ones_like(RRHO))
    FSRC = 'meff_span files (span_RANS50.npz, span_C193b.npz): sqrtk_mean, Leff_mean, cos_mean, rho_mean, regridded as the mu maps'
FFU, RKu, RLu, RCu, RMUu = FF, RK.copy(), RLn.copy(), RC.copy(), RMU.copy()      # unconditional point-file span means (previous maps)
ZSCF = np.interp(xf, xc[oks], zss[oks])
avail = avail & (ZF > 0)                             # wall nodes excluded: the first row is the first cell above the wall
COND = None
if MASK == '3d':
    xf, zf, Fc, nfc, nvc, CC = FC.get(S0)
    XF, ZF = np.meshgrid(xf, zf)
    with np.errstate(all='ignore'):
        RK = CC['k12R'] / CC['k12D']
        RLn = CC['LR'] / CC['LD']
        RC = 1.0 / CC['cD']
        RRHO = CC['rhoR'] / CC['rhoD']
        RMU = CC['muR'] / CC['muD']
        PROD = RRHO * RK * RLn * RC
    avail = (nfc >= 119) & (ZF > 0)
    GREYM = avail & (Fc < FMIN)
    FF = dict(cD=CC['cD'])
    COND = dict(f=Fc, n=nfc, nv=nvc, C=CC)
    FSRC = ('point file, strain mask in 3D before the span mean: fields interpolated onto the (x, y, z) grid as span_compute.py; mask = '
            f'interpolated |S*| D/U_inf > {S0} (193M, both rungs); conditional span means over the valid points; cells drawn where the valid '
            f'fraction f >= {FMIN}; cached span_cond_050.npz (_meff_factor_cond.py)')
else:
    GREYM = avail & (ZF > ZSCF[None, :])
print('FACTOR SOURCE:', FSRC)
INSIDE = (ZF < ZSCF[None, :]) if MASK == '2d' else (Fc >= FMIN)
below = avail & (XF >= XMIN) & (ZF > 0.0) & INSIDE & (RK > 0) & (RLn > 0) & (RC > 0) & (RMU > 0) & (RRHO > 0)
dln = np.log(PROD[below] / RMU[below])
print(f'product check (x/D 3-20, 0 < z < z_s, all ratios > 0, {below.sum()} cells): ln(product of span-mean factor ratios / span-mean mu ratio): '
      f'median {np.median(dln):+.3f}, median |.| {np.median(np.abs(dln)):.3f}, 90th percentile of |.| {np.percentile(np.abs(dln), 90):.3f}')
PCHK = dict(n=int(below.sum()), med=float(np.median(dln)), medabs=float(np.median(np.abs(dln))), p90=float(np.percentile(np.abs(dln), 90)))
# band where the span-mean mu ratio is within 10 % of 1
BAND_X = [5, 10, 15, 20]
BANDSTAT = {}
for xq in BAND_X:
    j = int(np.abs(xf - xq).argmin())
    sel_ = below[:, j] & (RMU[:, j] >= 0.9) & (RMU[:, j] <= 1.1)
    if sel_.any():
        BANDSTAT[xq] = dict(z0=zf[sel_].min(), z1=zf[sel_].max(), n=int(sel_.sum()), k=RK[sel_, j].mean(), L=RLn[sel_, j].mean(),
                            ic=RC[sel_, j].mean(), rho=RRHO[sel_, j].mean(), mu=RMU[sel_, j].mean())
        b_ = BANDSTAT[xq]
        print(f'band mu ratio 0.9-1.1 at x/D {xq}: z/D {b_["z0"]:.2f}-{b_["z1"]:.2f} ({b_["n"]} cells): mean ratios k^1/2 {b_["k"]:.2f}, L {b_["L"]:.2f}, '
              f'1/cos {b_["ic"]:.2f}, rho {b_["rho"]:.3f}, mu {b_["mu"]:.2f}')
ZPROF = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
PROF = {}
for xq in BAND_X:
    j = int(np.abs(xf - xq).argmin())
    rows = []
    for zq in ZPROF:
        i = int(np.abs(zf - zq).argmin())
        rows.append((zq, RK[i, j], RLn[i, j], RC[i, j], RRHO[i, j], PROD[i, j], RMU[i, j], FF['cD'][i, j]))
    PROF[xq] = rows
for r_ in PROF[10]:
    print('x/D 10 z/D %.1f: k^1/2 %.2f  L %.2f  1/cos %.2f  rho %.3f  product %.2f  mu ratio %.2f  (cos193 %.2f)' % r_)

ZX_HW = 0.1


def zline(R_, xg=None, zg=None, hw=ZX_HW):
    """z_x(x): height where, going up from the wall (z > 0), the ratio R_ (columns = x) reaches 1 for the first time (linear interpolation in z);
    light running median over +-hw in x/D (hw = 0.1 D: 5 columns of 0.05 D, only to remove cell-scale steps). Returns (raw, smoothed)."""
    xg = xf if xg is None else xg
    zg = zf if zg is None else zg
    raw = np.array([crossing(R_[:, j], np.ones(len(zg)), zg) for j in range(len(xg))])
    sm = np.full(len(xg), np.nan)
    for j in range(len(xg)):
        w = (np.abs(xg - xg[j]) <= hw + 1e-9) & np.isfinite(raw)
        if w.sum() >= 3:
            sm[j] = np.median(raw[w])
    return raw, sm


ANA = {}
REGG = (nfc >= 119) & (ZF > 0) & (ZF <= ZMAX + 1e-9) & (XF >= XMIN) & (XF <= 20) if MASK == '3d' else None
if MASK == '3d':
    # (1) sensitivity to f_min: grey area fraction and height of the upper edge (lowest z with f < f_min)
    ANA['fmin'] = {}
    patch = (XF >= 3) & (XF <= 5) & (ZF >= 1.0) & (ZF <= 1.8)           # top-left red patch of the earlier maps
    for fm in (0.02, 0.1, 0.2, 0.25, 0.3, 0.4, 0.42, 0.45, 0.5):
        gfrac = float(((Fc < fm) & REGG).sum() / REGG.sum())
        edges = []
        for xq in BAND_X:
            j = int(np.abs(xf - xq).argmin())
            idx = np.where((Fc[1:, j] < fm))[0]
            edges.append(float(zf[1:][idx[0]]) if len(idx) else np.nan)
        with np.errstate(all='ignore'):
            sat = {nm: int((avail & (Fc >= fm) & patch & (rr >= 8)).sum()) for nm, rr in (('L', RLn), ('1/cos', RC), ('mu', RMU))}
        ANA['fmin'][fm] = dict(grey=gfrac, edges=edges, sat=sat)
        print(f'f_min {fm}: grey fraction of the plotted area {100 * gfrac:.1f} %, upper edge z/D at x/D 5, 10, 15, 20: ' + ', '.join(f'{v:.2f}' for v in edges) +
              f' | cells with ratio >= 8 in the top-left patch (x/D 3-5, z/D 1-1.8): {sat}')
    # (2) conditional against unconditional maps below the old z_s (where both are defined)
    both = avail & (XF >= XMIN) & (ZF < ZSCF[None, :]) & (Fc >= FMIN)
    ANA['diff'] = {}
    for nm, new, old in (('mu', RMU, RMUu), ('k^1/2', RK, RKu), ('L', RLn, RLu), ('1/cos', RC, RCu)):
        with np.errstate(all='ignore'):
            m_ = both & (new > 0) & (old > 0)
            dd = np.abs(np.log(new[m_] / old[m_]))
        ANA['diff'][nm] = (int(m_.sum()), float(np.median(dd)), float(np.percentile(dd, 90)))
        print(f'conditional against unconditional ratio {nm}: |ln difference| median {np.median(dd):.3f}, 90th percentile {np.percentile(dd, 90):.3f} ({m_.sum()} cells below the old z_s)')
    # (3) heights of the ratio = 1 lines, new and old
    ANA['lines'] = {}
    for nm, new, old in (('mu', RMU, RMUu), ('k^1/2', RK, RKu)):
        zn, zo = zline(np.where(avail & (Fc >= FMIN), new, np.nan))[1], zline(np.where(nfc >= 119, old, np.nan))[1]
        ANA['lines'][nm] = ([float(np.interp(xq, xf, zn)) for xq in BAND_X], [float(np.interp(xq, xf, zo)) for xq in BAND_X])
        print(f'height of the {nm} ratio = 1 line at x/D 5, 10, 15, 20: new ' + ', '.join(f'{v:.3f}' for v in ANA['lines'][nm][0]) +
              ' | old ' + ', '.join(f'{v:.3f}' for v in ANA['lines'][nm][1]))
    ZXN_raw, ZXN = zline(np.where(avail & (Fc >= FMIN), RMU, np.nan))               # point file, conditional span means
    print(f'z_x (point file): max |raw - smoothed| = {np.nanmax(np.abs(ZXN_raw - ZXN)):.3f} D (window +-{ZX_HW} D)')
    # (4) mask of the bars against the mask of the maps; volume mean against the f-weighted mean of the conditional span means
    ANA['chk'] = {}
    for rk, (lo, hi) in (('near film', (3, 10)), ('far film', (10, 20.0001))):
        reg_ = REGG & (XF >= lo) & (XF < hi)
        nvg, nfg = int(nvc[reg_].sum()), int(nfc[reg_].sum())
        mk_ = MK[(rk, 'all')]
        fw = {}
        for key in ('muR', 'muD'):
            a_ = CC[key]
            m2 = reg_ & np.isfinite(a_)
            fw[key] = (float((Fc[m2] * a_[m2]).sum() / Fc[m2].sum()), float((nvc[m2] * a_[m2]).sum() / nvc[m2].sum()))
        ANA['chk'][rk] = dict(nvg=nvg, nfg=nfg, gridfrac=nvg / nfg, bars_n=mk_['n'], bars_frac=mk_['kept'], fw=fw, bars=(mk_['muRm'], mk_['muDm']))
        print(f'{rk}: grid valid span points {nvg} of {nfg} ({100 * nvg / nfg:.1f} %); bars kept points {mk_["n"]} (volume fraction {100 * mk_["kept"]:.1f} %); '
              f'f-weighted mean of the conditional means mu_t {fw["muR"][0]:.5f}, mu_eff {fw["muD"][0]:.5f} | direct masked mean (valid-count weighted) {fw["muR"][1]:.5f}, {fw["muD"][1]:.5f} '
              f'| bars (volume mean) {mk_["muRm"]:.5f}, {mk_["muDm"]:.5f}')

# ---------------------------------------------------------------- near-wall boundary z_x(x) and the continuous split
XQ = np.arange(XMIN, 20 + 1e-9, 0.25)
ZX_MS = None
XFP, ZFP = xf, zf                                    # point-file grid (kept for the comparison lines)
XZM = xc
if SRC == 'span':
    # full-pitch conditional span means of the meff_span files, both runs on the 193M grid (RANS interpolated, see _meff_span_split.py)
    SG = SS2.load()
    with np.errstate(all='ignore'):
        _rat = SG['R']['mu'] / SG['D']['mu']
        _fv = SG['f']
        _okc = np.isfinite(_rat) & np.isfinite(_fv) & (_fv >= FMIN)
        _, ZX_MS = zline(np.where(_okc, _rat, np.nan), SG['xc'], SG['zc'])
        RK, RLn = SG['R']['k'] / SG['D']['k'], SG['R']['L'] / SG['D']['L']
        RC, RRHO, RMU = 1.0 / SG['D']['c'], SG['R']['rho'] / SG['D']['rho'], SG['R']['mu'] / SG['D']['mu']
        PROD = RRHO * RK * RLn * RC
    xf, zf, Fc = SG['xc'], SG['zc'], SG['f']
    XZM = xf
    XF, ZF = np.meshgrid(xf, zf)
    nfc = np.where(np.isfinite(Fc), 999, 0)
    avail = np.isfinite(Fc) & (ZF > 0)
    GREYM = avail & (Fc < FMIN)
    AREAW = SG['area']
    FSRC = f'meff_span conditional span means (c_*_mean), full pitch, RANS interpolated onto the 193M grid; valid fraction f >= {FMIN}'
    print('FACTOR SOURCE:', FSRC)
if MASK == '3d':
    if SRC == 'span':
        ZLINE, ZSRC = (xf, ZX_MS), ('meff_span conditional span means (span_RANS50.npz, span_C193b.npz: c_mu_mean, c_frac_valid), full pitch, '
                                    'RANS bilinearly interpolated onto the 193M (x, z) bin centres')
    elif ZX_MS is not None:
        ZLINE, ZSRC = (xc, ZX_MS), 'meff_span conditional span means (span_RANS50.npz, span_C193b.npz: c_mu_mean, c_frac_valid), full pitch'
    else:
        ZLINE, ZSRC = (xf, ZXN), 'point file (pts.npz), conditional span means from _meff_factor_cond.py (meff_span conditional arrays not yet present)'
else:
    ZLINE, ZSRC = (xc[okx], zxs[okx]), 'unconditional span means (MASK = 2d)'
print('z_x SOURCE:', ZSRC)
_okz = np.isfinite(ZLINE[1])
ZXTAB = (XQ, np.interp(XQ, ZLINE[0][_okz], ZLINE[1][_okz]))
ZXCMP = [3, 4, 5, 7, 10, 15, 20]
print('z_x at x/D ' + ', '.join(str(v) for v in ZXCMP) + ': ' + ', '.join(f'{np.interp(v, *ZXTAB):.3f}' for v in ZXCMP))
if MASK == '3d':
    print('   point file     : ' + ', '.join(f'{np.interp(v, XFP[np.isfinite(ZXN)], ZXN[np.isfinite(ZXN)]):.3f}' for v in ZXCMP))
    if ZX_MS is not None:
        print('   meff_span      : ' + ', '.join(f'{np.interp(v, XZM[np.isfinite(ZX_MS)], ZX_MS[np.isfinite(ZX_MS)]):.3f}' for v in ZXCMP))
    CEN = np.arange(XMIN + 0.5, 19.5 + 1e-9, 0.25)               # slabs of 1 D (x/D XMIN-20), centres every 0.25 D
    TABC = [5, 7, 10, 15, 20]
    # primary split: point file, products exact per point, 3D strain mask (SP); cross-check: bin means of the meff_span files (SS2)
    cc_ = lambda c_, w_=1.0, **kw: SP.compute_curve(ZXTAB, c_, w_, **kw)  # noqa: E731
    ccs_ = lambda c_, w_=1.0, **kw: SS2.compute_curve(SG, ZXTAB, c_, w_, fmin=FMIN, **kw) if SRC == 'span' else cc_(c_, w_, **kw)  # noqa: E731
    CUR = cc_(CEN)
    CURS = cc_(TABC)
    CURA = cc_([5, 10, 15, 20], layer='above')
    FILMNW = cc_([(XMIN + 20) / 2], 20 - XMIN)[0]
    FILMA = cc_([(XMIN + 20) / 2], 20 - XMIN, layer='above')[0]
    XCUR = ccs_(CEN)
    XCURS = ccs_(TABC)
    XFILM = ccs_([(XMIN + 20) / 2], 20 - XMIN)[0]
    XFILMA = ccs_([(XMIN + 20) / 2], 20 - XMIN, layer='above')[0]
    print('continuous split: max |sum of terms - total| over the %d slabs: %.1e' % (len(CUR), max(abs(r['res']) for r in CUR)))
    for r in CURS:
        print('  x/D %4.1f near wall: total %+.3f terms ' % (r['x'], r['tot']) + ' '.join(f'{v:+.3f}' for v in r['p4']) +
              ' kept %.3f fneg %.4f' % (r['kept'], r['fneg']))

# hole opening from the diffuser/exit points near the wall
pd_ = np.load(M.SCR + 'pts.npz')
Xp = np.load(M.SCR + 'xyz.npz')['xyz']
names = [str(s) for s in pd_['regnames']]
mm = np.isin(pd_['reg'], [names.index('hole_diffuser'), names.index('exit')]) & (Xp[:, 2] > -0.03) & (Xp[:, 2] < 0.0)
H0, H1 = Xp[mm, 0].min(), Xp[mm, 0].max()
print(f'hole opening at the wall: x/D {H0:.2f} to {H1:.2f}')

# ---------------------------------------------------------------- figure
set_mpl_params(10.0 / 22.0, grid_on=False)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
with np.errstate(all='ignore'):
    ok = (Rm > 0) & (Dm > 0)
    lr = np.where(ok, np.log10(Rm / Dm), np.nan)
cell = np.isfinite(Rm) & np.isfinite(Dm)
ze2 = ze[:np.searchsorted(ze, ZMAX + 1e-9) + 1]
if ze2[-1] < ZMAX - 1e-9:
    ze2 = np.r_[ze2, ZMAX]
nzp = len(ze2) - 1
if SRC == 'span':
    xef, zef = SG['xe'], SG['ze'].copy()             # true bin edges of the 193M grid
    zef[0] = 0.0
else:
    xef = np.r_[xf[0] - 0.5 * (xf[1] - xf[0]), 0.5 * (xf[1:] + xf[:-1]), xf[-1] + 0.5 * (xf[-1] - xf[-2])]
    zef = np.r_[0.0, 0.5 * (zf[1:] + zf[:-1]), zf[-1] + 0.005]
zef1 = zef[1:].copy()
zef1[0] = 0.0                                        # wall row (z = 0) dropped; the first cell above the wall starts at the wall
ze2d = ze2[1:].copy()
ze2d[0] = 0.0

PANELS = os.environ.get('F14_PANELS', 'd')       # 'd' (default): mu panel only, one column; 'abcd': the four panels (k^1/2, L, 1/cos theta, mu)
# 20M LES panel (resolution ladder, owner 2026-11-27): F14_LES=1 (default) adds panel (b) below the RANS panel (a); F14_LES=0 gives the previous
# one-panel figure. The point file has no LES counterpart, so panel (b) is formed from the meff_span files: ratio of the conditional
# Delta-y-weighted span means cy_mu_mean, span_C20bW.npz (20M LES, mu_eff = resolved stress + mean WALE stress) over span_C193b.npz, both
# interpolated bilinearly (bin centres) onto the (x, z) grid of panel (a). Panel (a) is not changed. Masks (valid fraction f >= FMIN,
# manual jet-edge mask, ratio <= 0) and the z_x line are those of panel (a).
LES = os.environ.get('F14_LES', '1') == '1' and PANELS != 'abcd' and MASK == '3d' and SRC == 'pointfile'
RLES = None
if LES:
    from scipy.interpolate import RegularGridInterpolator as _RGI
    _fl, _fd, _fr = (np.load(SPAN_DIR + f'span_{c_}.npz') for c_ in ('C20bW', 'C193b', 'RANS50'))

    def _on_grid(f_, key='cy_mu_mean'):
        ip_ = _RGI((f_['z_D'], f_['x_D']), f_[key].T.astype(np.float64), bounds_error=False, fill_value=None, method='linear')
        return ip_(np.stack(np.meshgrid(np.clip(zf, f_['z_D'][0], f_['z_D'][-1]), np.clip(xf, f_['x_D'][0], f_['x_D'][-1]), indexing='ij'), -1))
    # Ratio formed on common 0.1 D x bins (owner 2026-11-27), as the table route: LES and 193M cy_mu_mean rebinned in x to 0.1 D by count-weighted
    # averages (weight = c_mu_N times the overlap fraction; the cy weight sum is not stored; the merged thin 193M bin at x/D 10.15 is replaced by
    # the interpolation of its neighbours), the LES interpolated in z onto the 193M levels, ratio formed on these bins, then mapped onto the grid.
    import _meff_span_split as _SS
    _k0 = dict(_SS.KEYS)
    _SS.KEYS = dict(mu='cy_mu_mean', rho='cy_rho_mean', k='cy_sqrtk_mean', L='cy_Leff_mean', c='cy_cos_mean')
    try:
        _G = _SS.load(case='C20bW')            # both runs on common 0.1 D x bins (count-weighted), LES interpolated in z onto the 193M levels
    finally:
        _SS.KEYS = _k0
    with np.errstate(all='ignore'):
        _rl = _G['R']['mu'] / _G['D']['mu']                                           # ratio on the 0.1 D bins [nz193, nx]
    _fl = dict(z_D=_G['zc'], x_D=_G['xc'])

    def _to_grid(a_, f_):
        ip_ = _RGI((f_['z_D'], f_['x_D']), a_, bounds_error=False, fill_value=None, method='linear')
        return ip_(np.stack(np.meshgrid(np.clip(zf, f_['z_D'][0], f_['z_D'][-1]), np.clip(xf, f_['x_D'][0], f_['x_D'][-1]), indexing='ij'), -1))
    with np.errstate(all='ignore'):
        _dm = _on_grid(_fd)
        RLES = _to_grid(_rl, _fl)
        RSPAN_R = _on_grid(_fr) / _dm                       # RANS by the previous span route (comparison only, not drawn)
        _cmp = avail & (Fc >= FMIN) & (XF >= XMIN) & (ZF <= ZMAX) & np.isfinite(RSPAN_R) & (RSPAN_R > 0) & (RMU > 0)
        _dl = np.abs(np.log(RSPAN_R[_cmp] / RMU[_cmp]))
    print(f'LES panel: cy_mu_mean ratio span_C20bW / span_C193b on the point-file grid; RANS by the same span route against the drawn RANS '
          f'(point file) over {_cmp.sum()} drawn cells: |ln ratio| median {np.median(_dl):.3f}, 90th percentile {np.percentile(_dl, 90):.3f}')
if PANELS == 'abcd':
    fig = plt.figure(figsize=(5.625, 5.0))          # printed size: figure* at 0.75 \textwidth (5.625 in), scale 1
    gs = fig.add_gridspec(4, 2, width_ratios=[1, 0.03], hspace=0.14, wspace=0.03, left=0.10, right=0.865, top=0.965, bottom=0.095)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[1, 0], sharex=axa)
    axc = fig.add_subplot(gs[2, 0], sharex=axa)
    axd = fig.add_subplot(gs[3, 0], sharex=axa)
    cax = fig.add_subplot(gs[0:4, 1])
else:
    # Flow-field map at full text width (owner 2026-11-27): 7.5 in (figure*, width=\textwidth, scale 1), 10 pt text.
    # F14_ZX = z stretch of the axes (1 = true aspect, equal scales in x and z; at most 2). Default: the version that reads better.
    ZX_ = float(os.environ.get('F14_ZX', '2'))
    AXL_, AXW_, GAP_, CBW_, AXB_, AXT_ = 0.55, 6.3, 0.12, 0.11, 0.45, 0.08
    AXH_ = AXW_ * ZMAX / (20.0 - XMIN) * ZX_
    FW_, FH_ = 7.5, AXB_ + AXH_ + AXT_
    GAPV_ = 0.12                                              # vertical gap between the RANS and LES panels
    if LES:
        FH_ = AXB_ + 2 * AXH_ + GAPV_ + AXT_
    fig = plt.figure(figsize=(FW_, FH_))
    if LES:
        axd = fig.add_axes([AXL_ / FW_, (AXB_ + AXH_ + GAPV_) / FH_, AXW_ / FW_, AXH_ / FH_])      # (a) RANS 50M
        axl = fig.add_axes([AXL_ / FW_, AXB_ / FH_, AXW_ / FW_, AXH_ / FH_], sharex=axd)            # (b) 20M LES
        cax = fig.add_axes([(AXL_ + AXW_ + GAP_) / FW_, AXB_ / FH_, CBW_ / FW_, (2 * AXH_ + GAPV_) / FH_])
    else:
        axd = fig.add_axes([AXL_ / FW_, AXB_ / FH_, AXW_ / FW_, AXH_ / FH_])
        cax = fig.add_axes([(AXL_ + AXW_ + GAP_) / FW_, AXB_ / FH_, CBW_ / FW_, AXH_ / FH_])
    print(f'figsize {FW_} x {FH_:.2f} in, axes {AXW_} x {AXH_:.2f} in each ({2 if LES else 1} panel(s)), z stretch {ZX_:g}')


def wall(ax):
    ax.plot([XMIN, 20], [0, 0], 'k', lw=1.2, clip_on=False)


def frame(ax):
    wall(ax)
    ax.set_xlim(XMIN, 20)
    ax.set_ylim(0, ZMAX)
    ax.set_ylabel(r'$z/D$')
    ax.plot(ZLINE[0], ZLINE[1], 'k', lw=0.6, zorder=5)                       # z_x(x): one thin line on every map, no label


GREYC = matplotlib.colors.ListedColormap(['#d9d9d9'])
tag_kw = dict(ha='right', va='center', bbox=dict(fc='w', ec='none', pad=1.0, alpha=0.8))


def grey(ax, m, xedge, zedge, n_z):
    ax.pcolormesh(xedge, zedge, np.ma.masked_where(~m, np.ones(m.shape)), cmap=GREYC, rasterized=True)


def factor_map(ax, ratio, label):
    with np.errstate(all='ignore'):
        defined = np.isfinite(ratio) & (ratio > 0)
        lg = np.where(defined, np.log10(ratio), np.nan)
    zline_ = JET_MASK_Z0 + (JET_MASK_Z1 - JET_MASK_Z0) * (XF - JET_MASK_X0) / (JET_MASK_X1 - JET_MASK_X0)
    JETM = (XF >= JET_MASK_X0 - 1e-9) & (XF <= JET_MASK_X1) & (ZF > zline_)                  # manual jet-edge mask, drawing only
    ab = GREYM | (JETM & avail)                                     # no strain: above z_s (2d) or valid fraction f < f_min (3d); plus the jet-edge mask
    greyc = avail & (ab | ~defined)
    lg = np.where(ab, np.nan, lg)
    grey(ax, greyc[1:], xef, zef1, len(zf))
    im_ = ax.pcolormesh(xef, zef1, np.ma.masked_invalid(np.where(avail, lg, np.nan))[1:], cmap='bwr', norm=Normalize(-1, 1), rasterized=True)
    frame(ax)
    ax.text(0.992, 0.08, label, transform=ax.transAxes, ha='right', va='bottom',
            bbox=dict(fc='w', ec='none', pad=1.0, alpha=0.75))
    return im_, greyc


lab = dict(a=r'$k^{1/2}$', b=r'$L$', c=r'$1/\cos\theta$', d=r'$\mu$')
DRAWMED = SRC == 'span' and os.environ.get('F14_DRAWMED', '0') == '1'   # 3-column running median in x on the drawn ratio maps only (striping x/D 8-10.2)


def med3(a):
    with np.errstate(all='ignore'):
        st = np.stack([np.roll(a, 1, axis=1), a, np.roll(a, -1, axis=1)])
        out = np.nanmedian(st, axis=0)
    out[:, [0, -1]] = a[:, [0, -1]]
    return np.where(np.isfinite(a), out, np.nan)


if DRAWMED:
    _zraw = zline(np.where(np.isfinite(RMU) & (Fc >= FMIN), RMU, np.nan), xf, zf)[1]
    _zmed = zline(np.where(np.isfinite(RMU) & (Fc >= FMIN), med3(RMU), np.nan), xf, zf)[1]
    _m = (xf >= 3.5) & np.isfinite(_zraw) & np.isfinite(_zmed)
    print(f'3-column running median in x of the ratio maps: z_x(median-filtered) - z_x max |diff| = {np.abs(_zraw - _zmed)[_m].max():.4f} D (x/D 3.5-20)')
    RK, RLn, RC, RMU = med3(RK), med3(RLn), med3(RC), med3(RMU)          # drawn maps only; z_x and the split use the unfiltered ratios
if PANELS == 'abcd':
    im, g_a = factor_map(axa, RK, lab['a'])
    factor_map(axb, RLn, lab['b'])
    factor_map(axc, RC, lab['c'])
if MASK == '3d':
    _imd = factor_map(axd, RMU, '(a) RANS 50M' if LES else lab['d'])
    if PANELS != 'abcd':
        im = _imd[0]
    if LES:
        if os.environ.get('F14_LES_MED', '0') == '1':        # 3-column running median in x on panel (b) only (off by default since the 0.1 D rebin; the 0.05 D LES bins hold one mesh x line each and alternate by about 10 %)
            RLES = med3(RLES)
        factor_map(axl, RLES, '(b) 20M LES')
        plt.setp(axd.get_xticklabels(), visible=False)
        axl.set_xlabel(r'$x/D$', labelpad=1)
    zhalf = [float(zf[1:][np.where(Fc[1:, int(np.abs(xf - xq).argmin())] < 0.5)[0][0]]) if (Fc[1:, int(np.abs(xf - xq).argmin())] < 0.5).any() else np.nan
             for xq in BAND_X]
else:
    # (d) ratio of the span-mean mu (meff_span source)
    lrd = np.where(LOW, np.nan, lr)
    axd.pcolormesh(xe, ze2d, np.ma.masked_where(~(cell & (LOW | ~ok))[1:nzp], np.ones((nzp - 1, nx))), cmap=GREYC, rasterized=True)
    axd.pcolormesh(xe, ze2d, np.ma.masked_invalid(lrd[1:nzp]), cmap='bwr', norm=Normalize(-1, 1), rasterized=True)
    frame(axd)
    axd.text(0.992, 0.08, lab['d'], transform=axd.transAxes, ha='right', va='bottom', bbox=dict(fc='w', ec='none', pad=1.0, alpha=0.75))
if not LES:
    axd.set_xlabel(r'$x/D$', labelpad=1)
if PANELS == 'abcd':
    for ax in (axa, axb, axc):
        plt.setp(ax.get_xticklabels(), visible=False)
cb = fig.colorbar(im, cax=cax)
rt = [0.1, 0.2, 0.5, 1, 2, 5, 10]
_short = PANELS != 'abcd' and AXH_ < 1.2                  # true-aspect single panel: fewer ticks, short label
if _short:
    rt = [0.1, 1, 10]
cb.set_ticks(np.log10(rt))
cb.set_ticklabels([str(v) for v in rt])
cb.set_label('ratio RANS / 193M of the span means' if PANELS == 'abcd' else ('ratio' if _short else ('ratio to 193M' if LES else 'RANS / 193M')), labelpad=2)
# panel letters only (the caption explains the grey); the single-panel figure has none
if PANELS == 'abcd':
    for ax, t in ((axa, '(a)'), (axb, '(b)'), (axc, '(c)'), (axd, '(d)')):
        ax.text(-0.095, 1.0, t, transform=ax.transAxes, ha='left', va='bottom')
fig.savefig(OUT, dpi=600)
print('saved', OUT)
_gx = (nfc >= 119) & (ZF > 0) & (ZF <= ZMAX + 1e-9) & (XF >= XMIN) & (XF <= 20)
print(f'grey fraction of the plotted area (x/D {XMIN:g}-20, f_min {FMIN}): {100 * ((Fc < FMIN) & _gx).sum() / _gx.sum():.1f} % of bins')
_zl = JET_MASK_Z0 + (JET_MASK_Z1 - JET_MASK_Z0) * (XF - JET_MASK_X0) / (JET_MASK_X1 - JET_MASK_X0)
_jm = (XF >= JET_MASK_X0 - 1e-9) & (XF <= JET_MASK_X1) & (ZF > _zl) & _gx
_vis = _gx & (Fc >= FMIN) & ~_jm
_nsat = int(((_gx & (Fc >= FMIN)) & (RMU >= 3)).sum())
_nsat_h = int((_jm & (Fc >= FMIN) & (RMU >= 3)).sum())
print(f'without the manual mask: drawn cells with mu ratio >= 3 in the frame (x/D {XMIN:g}-20): {_nsat}, of them inside the mask: {_nsat_h}; '
      f'max drawn ratio inside the mask {np.nanmax(np.where(_jm & (Fc >= FMIN), RMU, np.nan)):.2f}')
print(f'manual jet-edge mask: line ({JET_MASK_X0}, {JET_MASK_Z0}) to ({JET_MASK_X1}, {JET_MASK_Z1}); cells hidden that were drawn: {(_jm & (Fc >= FMIN)).sum()}; '
      f'grey fraction with the mask: {100 * ((Fc < FMIN) | _jm)[_gx].sum() / _gx.sum():.1f} % of bins')
for _nm, _r in (('k^1/2', RK), ('L', RLn), ('1/cos', RC), ('mu', RMU)):
    print(f'  remaining drawn cells with ratio >= 3 in x/D 3-6, z/D 0.9-2, panel {_nm}: {int((_vis & (XF <= 6) & (ZF > 0.9) & (_r >= 3)).sum())}')
if SRC == 'span':
    _gm = (avail & (ZF <= ZMAX + 1e-9) & (XF >= XMIN) & (XF <= 20))
    _gr = _gm & (GREYM | ~(np.isfinite(RMU) & (RMU > 0)))
    print(f'   area-weighted (bin area dx dz), all grey cells (f < f_min or ratio undefined): {100 * AREAW[_gr].sum() / AREAW[_gm].sum():.1f} %; '
          f'f < f_min only: {100 * AREAW[_gm & GREYM].sum() / AREAW[_gm].sum():.1f} %; bins {_gm.sum()}')
if ALT:
    sys.exit(0)                                       # alternative version: figure only

# ---------------------------------------------------------------- review sheet: raw factor fields
if SRC == 'pointfile':                              # raw point-file factor sheet (review only; not drawn for the meff_span source)
    FF = FFU                                        # review sheet: unconditional raw span means
    sel_ = avail
    figr, axr = plt.subplots(3, 2, figsize=(7.0, 6.6), sharex=True, gridspec_kw=dict(hspace=0.25, wspace=0.12, left=0.07, right=0.88, top=0.98,
                                                                                   bottom=0.07))
    rawf = [(FF['k12R'] / UINF, r'$k^{1/2}/U_\infty$, RANS'), (FF['k12D'] / UINF, r'$k^{1/2}/U_\infty$, 193M'),
            (FF['LR'] / DD, r'$L_R/D$, RANS $\mu_t/(\rho k^{1/2})$'), (FF['LD'] / DD, r'$L_\mathrm{eff}/D$, 193M')]
    for (axr_, (arr, lab)) in zip(axr.ravel()[:4], rawf):
        a_ = np.where(sel_ & (arr > 0), arr, np.nan)
        if axr_ in (axr[0, 0], axr[0, 1]):
            vm = (1e-3, 0.3)
        else:
            vm = (1e-3, 0.5)
        ii = axr_.pcolormesh(xef, zef, np.ma.masked_invalid(a_), cmap='inferno', norm=LogNorm(*vm), rasterized=True)
        axr_.text(19.8, 1.85, lab, ha='right', va='center', fontsize=7, color='w')
        figr.colorbar(ii, ax=axr_, pad=0.01, aspect=15)
    a_ = np.where(sel_, FF['cD'], np.nan)
    ii = axr[2, 0].pcolormesh(xef, zef, np.ma.masked_invalid(a_), cmap='bwr', norm=Normalize(-1, 1), rasterized=True)
    axr[2, 0].text(19.8, 1.85, r'$\cos\theta$, 193M', ha='right', va='center', fontsize=7)
    figr.colorbar(ii, ax=axr[2, 0], pad=0.01, aspect=15)
    a_ = np.where(sel_ & (FF['rhoR'] > 0), FF['rhoR'] / FF['rhoD'], np.nan)
    ii = axr[2, 1].pcolormesh(xef, zef, np.ma.masked_invalid(np.log10(a_)), cmap='bwr', norm=Normalize(-1, 1), rasterized=True)
    axr[2, 1].text(19.8, 1.85, r'$\log_{10}(\rho_\mathrm{RANS}/\rho_{193\mathrm{M}})$', ha='right', va='center', fontsize=7)
    figr.colorbar(ii, ax=axr[2, 1], pad=0.01, aspect=15)
    for axr_ in axr.ravel():
        axr_.set_xlim(2.8, 20)
        axr_.set_ylim(0, 2)
    for axr_ in axr[:, 0]:
        axr_.set_ylabel(r'$z/D$')
    for axr_ in axr[2]:
        axr_.set_xlabel(r'$x/D$')
    figr.savefig(RAW, dpi=200)
    print('saved', RAW)

# ---------------------------------------------------------------- table
L = []
P = L.append
P('# F14 (fig_meff_maps): factor maps and exact factor split, RANS 50M SST against 193M (b010)')
P(f'BOUNDARY = {BOUNDARY}, WEIGHTING = {WEIGHTING}')
P('Source of the span means of mu (map d): ' + S['src'])
P('Source of the span means of the factors (maps a-c): point file (pts.npz), interpolated as span_compute.py, plain mean over the points '
  'of the span line, complete span only (>= 119 of 121 points, x/D >= 2.8); cached span_factors.npz (_meff_factor_span.py).')
P('  Arrays to request from the other session (per rung, same bins as the mu files, plain span mean over all grid points of the pitch, '
  'wall nodes included, same N / N_nonfinite): rho, k^(1/2) (mean of the point-wise square root of k), L_eff (193M: k^(1/2)|b|/|S*|; '
  'RANS: mu_t/(rho k^(1/2))), cos(theta) (193M; RANS = 1), L_al = L_eff cos(theta) (193M), and k itself.')
P('  Lengths: 193M L_eff = L of the point file, mu_eff = rho k^(1/2) L_eff cos(theta); RANS L_R = mu_t/(rho_R k_R^(1/2)).')
P(f'Colour scale of the ratio maps: log10 ratio -1 to +1 (0.1 to 10); (c) is drawn as 1/cos(theta) (RANS cos = 1) so that the ratios multiply.')
P('')
P('z_x(x): lowest height where span-mean RANS mu_t reaches the 193M span-mean mu_eff (RANS lower below, higher above), z/D > 0, '
  'running median over +-0.25 x/D.')
P('  x/D        ' + ' '.join(f'{v:7.1f}' for v in XT0))
P('  z_x/D, all cells (old)       ' + ' '.join(f'{v:7.3f}' for v in ZXT_OLD))
P('  z_x/D, defined cells (new)   ' + ' '.join(f'{v:7.3f}' for v in ZXT))
P('  table used by the strain-weighted split before: ' + ' '.join(f'{v:.2f}' for v in ZT0))
P(f'  max difference {dmax:.3f} D; strain-weighted split recomputed with the new z_x: {"yes" if RERUN else "no"}.')
P('z_s(x): lowest height above the wall where the 193M span-mean |S*| D/U_inf falls below 0.05, running median over +-0.25 x/D.')
P('  z_s/D        ' + ' '.join(f'{v:7.3f}' for v in ZST))
P('')
P(f'MASK = {MASK}, FMIN = {FMIN}, LAYOUT = {LAYOUT}')
if MASK == '3d':
    P('Method (MASK = 3d): the strain mask acts point by point in 3D before the span average. All fields are interpolated onto the regular (x, y, z) grid')
    P('  (linear in (x, y) on the common z levels, linear in z, as span_compute.py); the mask is computed on that grid from the INTERPOLATED point-wise')
    P(f'  193M |S*| D/U_inf, valid = interpolated strain > {S0} (a threshold on the interpolated strain, no 0/1 flag is interpolated); same mask for both rungs.')
    P('  At each (x, z) the conditional span mean is the mean over the valid points of the span line (finite for the field); f = valid / finite-strain points.')
    P(f'  A cell is drawn where f >= {FMIN} (complete span, >= 119 of 121 points); grey otherwise ("no mean strain") and where a mean is <= 0.')
    P('Sensitivity to f_min (plotted area x/D 3-20, z/D 0-2, wall row excluded; upper edge = lowest z where f < f_min; nan = f stays >= f_min up to z/D 2):')
    P(f'  {"f_min":>6s} {"grey area":>10s} {"edge x/D 5":>11s} {"10":>6s} {"15":>6s} {"20":>6s}')
    for fm, r_ in ANA['fmin'].items():
        P(f'  {fm:6.2f} {100 * r_["grey"]:9.1f}% ' + ' '.join(f'{v:6.2f}' for v in r_['edges']))
    P('  The f = 0.5 line (half the span strained) is drawn dashed on map (d); its height at x/D 5, 10, 15, 20: ' + ', '.join(f'{v:.2f}' for v in zhalf))
    P('Conditional against unconditional ratio maps (point-file span means) over the cells below the old z_s (span-mean strain boundary) with f >= FMIN:')
    for nm, (n_, med_, p90_) in ANA['diff'].items():
        P(f'  {nm:6s} |ln(new/old)| median {med_:.3f}, 90th percentile {p90_:.3f} ({n_} cells)')
    P('Height of the ratio = 1 line (lowest upward crossing, running median +-0.25 in x/D) at x/D 5, 10, 15, 20:')
    for nm, (zn_, zo_) in ANA['lines'].items():
        P(f'  {nm:6s} new ' + ' '.join(f'{v:.3f}' for v in zn_) + '   old ' + ' '.join(f'{v:.3f}' for v in zo_))
    P('Bars (point-wise 3D mask on the 50M-grid points) against the maps (interpolated grid): same definition (|S*| D/U_inf > 0.05, 193M, both rungs).')
    P('  Counts differ because the sets are different samplings (points of the mesh against grid points of the span lines); volume means:')
    for rk_, c_ in ANA['chk'].items():
        P(f'  {rk_}: grid valid span points {c_["nvg"]} of {c_["nfg"]} ({100 * c_["gridfrac"]:.1f} %); bars kept points {c_["bars_n"]} (volume fraction {100 * c_["bars_frac"]:.1f} %)')
        P(f'    mean mu_t/(rho_inf U_inf D): f-weighted conditional span means {c_["fw"]["muR"][0]:.5f}, direct masked mean {c_["fw"]["muR"][1]:.5f}, bars (volume mean) {c_["bars"][0]:.5f}')
        P(f'    mean mu_eff/(rho_inf U_inf D): f-weighted {c_["fw"]["muD"][0]:.5f}, direct {c_["fw"]["muD"][1]:.5f}, bars {c_["bars"][1]:.5f}')
    P('  f-weighted and direct masked means agree to rounding (same grid); the bars use mesh volumes as weights, hence the small difference.')
P('Product check: ln(rho x k^(1/2) x L x 1/cos ratios of the span means / ratio of the span-mean mu), point-file span means, '
  f'x/D 3-20, 0 < z < z_s, cells with all ratios > 0 ({PCHK["n"]} cells): median {PCHK["med"]:+.3f}, median |.| {PCHK["medabs"]:.3f}, '
  f'90th percentile of |.| {PCHK["p90"]:.3f}.')
P('')
P('Band where the span-mean mu ratio RANS/193M lies in 0.9-1.1 (cells below z_s; mean of the ratios over the band cells in the column):')
P(f'  {"x/D":>5s} {"z range":>12s} {"cells":>6s} {"k^1/2":>7s} {"L":>7s} {"1/cos":>7s} {"rho":>7s} {"mu":>7s}')
for xq, b_ in BANDSTAT.items():
    P(f'  {xq:5d} {b_["z0"]:5.2f}-{b_["z1"]:<6.2f} {b_["n"]:6d} {b_["k"]:7.2f} {b_["L"]:7.2f} {b_["ic"]:7.2f} {b_["rho"]:7.3f} {b_["mu"]:7.2f}')
P('')
P('Factor ratios through the height (span means, RANS/193M, nearest cell; product = rho x k^1/2 x L x 1/cos):')
for xq in BAND_X:
    P(f'  x/D {xq}:  {"z/D":>5s} {"k^1/2":>7s} {"L":>7s} {"1/cos":>7s} {"rho":>7s} {"product":>8s} {"mu":>7s} {"cos193":>7s}')
    for r_ in PROF[xq]:
        P('         ' + f'{r_[0]:5.1f} {r_[1]:7.2f} {r_[2]:7.2f} {r_[3]:7.2f} {r_[4]:7.3f} {r_[5]:8.2f} {r_[6]:7.2f} {r_[7]:7.2f}')
P('')
P('=' * 100)
P('BARS, strain as a mask only (WEIGHTING = mask): cells with the point-wise 193M |S*| D/U_inf > 0.05 kept, plain volume weights, wall nodes excluded,')
P('z <= 2 D, x/D 3-10 near film, 10-20 far film. Primary: total = ln(sum V mu_t,RANS / sum V mu_eff,193M), exact Shapley split into rho, k^1/2, L_eff, cos')
P('(RANS: L_R = mu_t/(rho k^1/2), cos = 1). Residual = sum of terms minus total. Kept = volume fraction of the (wall-node-free) region or band kept by the mask;')
P('mean top = mean over x/D bins of 0.25 of the largest z of the kept points; fneg = volume share of the kept cells with mu_eff,193M <= 0.')
P(f'  {"region / band":22s} {"total":>7s} {"rho":>7s} {"k^1/2":>7s} {"L_eff":>7s} {"cos":>7s} {"residual":>9s} {"kept":>6s} {"mean top":>9s} {"fneg":>6s} {"cells":>8s}')
for k, r in MK.items():
    P(f'  {k[0] + " / " + k[1]:22s} {r["tot"]:+7.3f} ' + ' '.join(f'{v:+7.3f}' for v in r['p4']) +
      f' {r["res"]:+9.1e} {r["kept"]:6.3f} {r["top"]:9.2f} {r["fneg"]:6.3f} {r["n"]:8d}')
P('  (z/D < 0.6 and z/D >= 0.6 are sub-bands of the kept cells; their kept fraction refers to the volume of the band.)')
P('Secondary: volume mean of the point-wise ln ratios over the kept cells with mu_eff,193M > 0 (exactly additive term by term; excluded = fneg above),')
P('and ratios RANS/193M of the volume-weighted medians per factor (cos: 1/median cos):')
P(f'  {"region / band":22s} {"ln mean: rho":>12s} {"k^1/2":>7s} {"L":>7s} {"1/cos":>7s} {"sum":>7s} | {"median ratios: rho":>18s} {"k^1/2":>7s} {"L":>7s} {"1/cos":>7s} {"mu":>7s}')
for k, r in MK.items():
    lm, md = r['lnmean'], r['med']
    P(f'  {k[0] + " / " + k[1]:22s} {lm["rho"]:12.3f} {lm["k"]:7.3f} {lm["L"]:7.3f} {lm["c"]:7.3f} {lm["tot"]:7.3f} | {md["rho"]:18.3f} {md["k"]:7.3f} {md["L"]:7.3f} {md["c"]:7.3f} {md["mu"]:7.3f}')
P('')
P('=' * 100)
P('Previous definition for the record (strain-squared weighted regional least squares, five players, layers by BOUNDARY):')
P('z_x(x): lowest height where span-mean RANS mu_t reaches the 193M span-mean mu_eff.')
P('Band ratios: ratio of band means (bin-area weighted mean of the RANS span mean / mean of the 193M span mean over the cells of the band; '
  'wall row excluded), x/D 3-10 near film, 10-20 far film; defined cells only (span-mean |S*| D/U_inf >= 0.05); "all cells" = no strain exclusion.')
P(f'  {"band":28s} {"mean RANS":>10s} {"mean 193M":>10s} {"ratio":>7s} {"all cells":>10s} {"cells":>6s}')
for k in TKEYS:
    b = BANDS[k]
    P(f'  {k[0] + " / " + k[1]:28s} {b["mR"]:10.5f} {b["mD"]:10.5f} {b["ratio"]:7.2f} {b["old"]:10.2f} {b["n"]:6d}')
P('')
P('Grey area of the maps (cells with span-mean |S*| D/U_inf < 0.05; owner exception to the no-masking rule, as for the Pr_t map); '
  'fraction of the plotted area of each layer (z <= 2 D, wall row excluded):')
P(f'  {"sub-region":28s} {"low strain":>11s} {"+ mean<=0":>11s} {"alt. definition":>16s} {"criterion & alt.":>17s}')
for k in TKEYS:
    g = GREY[k]
    P(f'  {k[0] + " / " + k[1]:28s} {100 * g["low"]:10.1f}% {100 * g["merged"]:10.1f}% {100 * g["alt"]:15.1f}% {100 * g["both"]:16.1f}%')
P('  alt. definition: fraction of span points with |S*| D/U_inf > 0.05 below one half (interpolated indicator).')
P(f'  x/D 3-20 plotted area: criterion {100 * (AREA * LOW)[areaw].sum() / AREA[areaw].sum():.1f} %, alternative '
  f'{100 * (AREA * LOWALT)[areaw].sum() / AREA[areaw].sum():.1f} %, symmetric difference {100 * (AREA * (LOW ^ LOWALT))[areaw].sum() / AREA[areaw].sum():.1f} %.')
P('')
P('Exact split (strain-squared weighted, Shapley over five factors), terms in ln(RANS/193M), wall nodes excluded, z <= 2 D:')
P(f'  {"sub-region":28s} {"total":>7s} {"rho":>7s} {"k^1/2":>7s} {"L_eff":>7s} {"cos":>7s} {"S-wt":>7s} {"resid":>8s}')
for k in TKEYS:
    r = SPL[k]
    P(f'  {k[0] + " / " + k[1]:28s} {r["tot"]:+7.3f} ' + ' '.join(f'{v:+7.3f}' for v in r['p5']) + f' {r["res"]:+8.0e}')
P('')
P('Plain volume means (no weighting, no filter, same cells as the strain-weighted split), shares and zero-or-negative fraction:')
P(f'  {"sub-region":28s} {"ln plain":>9s} {"ratio":>7s} {"vol share":>10s} {"strain-wt share":>16s} {"vol frac mu_eff<=0":>19s} {"cells":>8s}')
for k in TKEYS:
    r = SPL[k]
    P(f'  {k[0] + " / " + k[1]:28s} {r["plain"]:+9.3f} {np.exp(r["plain"]):7.2f} {r["vshare"]:10.3f} {r["sshare"]:16.3f} {r["fneg"]:19.3f} {r["n"]:8d}')
P('')
P('Ratio of volume-weighted medians of mu (RANS mu_t / 193M mu_eff, point-wise values); kept = volume fraction kept by the filter |S*|D/U > 0.05:')
P(f'  {"sub-region":28s} {"unfiltered":>11s} {"|S*|D/U>0.05":>13s} {"kept":>6s}')
for k in TKEYS:
    r = SPL[k]
    P(f'  {k[0] + " / " + k[1]:28s} {r["med_all"]:11.3f} {r["med_S05"]:13.3f} {r["keep_S05"]:6.3f}')
P('')
P('=' * 100)
P('NEAR-WALL LAYER (owner decision): z_x(x) = height where, going up from the wall, the ratio of the conditional span-mean eddy viscosity')
P(f'RANS/193M first reaches 1 (3D strain mask, valid fraction f >= {FMIN}); source: {ZSRC}; running median over +-{ZX_HW} D to remove cell-scale steps.')
P('  z_x/D at x/D ' + ' '.join(f'{v:6d}' for v in ZXCMP))
P('    used (authority)  ' + ' '.join(f'{np.interp(v, *ZXTAB):6.3f}' for v in ZXCMP))
if MASK == '3d':
    P('    point file        ' + ' '.join(f'{np.interp(v, XFP[np.isfinite(ZXN)], ZXN[np.isfinite(ZXN)]):6.3f}' for v in ZXCMP))
    if ZX_MS is not None:
        P('    meff_span         ' + ' '.join(f'{np.interp(v, XZM[np.isfinite(ZX_MS)], ZX_MS[np.isfinite(ZX_MS)]):6.3f}' for v in ZXCMP))
    else:
        P('    meff_span         not available yet (no c_mu_mean in span_RANS50.npz / span_C193b.npz)')
    P('    meff_span, unconditional full-pitch means (no mask, for reference) ' + ' '.join(f'{np.interp(v, xc[okx], zxs[okx]):6.3f}' for v in ZXCMP))
    P(f'  f_min = {FMIN} chosen as the smallest value (steps of the grid: 0.25, 0.3, 0.4, 0.42, 0.45, 0.5) at which the saturated top-left patch is gone from maps (b)-(d):')
    P(f'  {"f_min":>6s} {"grey":>7s}  cells with ratio >= 8 in x/D 3-5, z/D 1-1.8 (L, 1/cos, mu)')
    for fm, r_ in ANA['fmin'].items():
        P(f'  {fm:6.2f} {100 * r_["grey"]:6.1f}%  {r_["sat"]}')
    P('Continuous split in the near-wall layer (wall nodes excluded, z < z_x), slabs of 1 D in x (|x - centre| < 0.5 D), centres every 0.25 D from 4.0 to 19.5 (x/D 3.5-20);')
    P('point-wise 3D mask |S*| D/U_inf > 0.05, plain volume statistics; the x/D 20 row is the half slab 19.5-20.')
    P(f'  max |sum of the four terms - total| over all slabs: {max(abs(r["res"]) for r in CUR):.1e}')
    P(f'  {"x/D":>5s} {"total":>7s} {"rho":>7s} {"k^1/2":>7s} {"L_eff":>7s} {"cos":>7s} {"kept":>6s} {"mu_eff<=0":>10s}  ratios of the volume means k^1/2, L, 1/cos')
    for r in CURS:
        vm = r['vmean']
        P(f'  {r["x"]:5.1f} {r["tot"]:+7.3f} ' + ' '.join(f'{v:+7.3f}' for v in r['p4']) + f' {r["kept"]:6.3f} {r["fneg"]:10.4f}  '
          f'{vm["k"][0] / vm["k"][1]:.2f}  {vm["L"][0] / vm["L"][1]:.2f}  {1.0 / vm["c"][1]:.2f}')
    r = FILMNW
    P(f'  film mean x/D 3.5-20: {r["tot"]:+7.3f} ' + ' '.join(f'{v:+7.3f}' for v in r['p4']) + f' {r["kept"]:6.3f} {r["fneg"]:10.4f}')
    P('  (rows above: point file, products exact per point, so the four terms close on the total; layer = wall to the full-pitch z_x of the table above.)')
    P('  film mean above z_x, x/D 3.5-20 (point file): mu ratio %.3f, k^1/2 %.3f, L_eff %.3f, 1/cos %.3f (gap %+.1e)' % (
        np.exp(FILMA['tot']), FILMA['vmean']['k'][0] / FILMA['vmean']['k'][1], FILMA['vmean']['L'][0] / FILMA['vmean']['L'][1],
        1 / FILMA['vmean']['c'][1], FILMA['res']))
    if SRC == 'span':
        P(f'CROSS-CHECK, bin means of the meff_span files (f_min {FMIN}): weights dx dz f per 193M bin; RANS interpolated bilinearly onto the 193M bin centres;')
        P('  total = ln(sum w mu_t / sum w mu_eff) of the conditional means; the four terms are the Shapley split of ln(sum w prod_RANS / sum w prod_193M),')
        P('  prod = rho k^1/2 L_eff cos of the bin means. The mean of a product is not the product of the means: gap = sum of the terms - total.')
        P(f'  {"x/D":>5s} {"total":>7s} {"rho":>7s} {"k^1/2":>7s} {"L_eff":>7s} {"cos":>7s} {"gap":>7s}   point-file total')
        for r, rp in zip(XCURS, CURS):
            P(f'  {r["x"]:5.1f} {r["tot"]:+7.3f} ' + ' '.join(f'{v:+7.3f}' for v in r['p4']) + f' {r["res"]:+7.3f}   {rp["tot"]:+7.3f}')
        r = XFILM
        P(f'  film mean x/D 3.5-20: {r["tot"]:+7.3f} ' + ' '.join(f'{v:+7.3f}' for v in r['p4']) + f' {r["res"]:+7.3f}')
        P(f'  max |gap| over the slabs: {max(abs(r["res"]) for r in XCUR):.3f}; slabs with |gap| > 0.02: {sum(abs(r["res"]) > 0.02 for r in XCUR)} of {len(XCUR)}')
        P('  film mean above z_x (bin means): mu ratio %.3f, k^1/2 %.3f, L_eff %.3f, 1/cos %.3f (gap %+.3f)' % (
            np.exp(XFILMA['tot']), XFILMA['vmean']['k'][0] / XFILMA['vmean']['k'][1], XFILMA['vmean']['L'][0] / XFILMA['vmean']['L'][1],
            1 / XFILMA['vmean']['c'][1], XFILMA['res']))
    P('Layer above z_x (z >= z_x up to the strain-mask edge, z <= 2 D), slabs of 1 D (x/D 20: half slab):')
    for r in CURA:
        P(f'  {r["x"]:5.1f} {r["tot"]:+7.3f} ' + ' '.join(f'{v:+7.3f}' for v in r['p4']) + f' {r["kept"]:6.3f} {r["fneg"]:10.4f}')
P('')
P('Source used: ' + S['short'])
open(TAB, 'w').write('\n'.join(L) + '\n')
print('saved', TAB)
# full curves for the flux script (x/D 3 to 20 every 0.25)
xq_ = XQ
DDIR = os.path.dirname(TAB)
for fn, (xx_, zz_), nm in (('near_wall_boundary.txt', (xc[oks], zss[oks]), 'z_s/D (|S*| D/U_inf = 0.05 contour of the 193M span-mean strain)'),):
    with open(os.path.join(DDIR, fn), 'w') as fh:
        fh.write(f'# x/D  {nm}; source: {S["short"]}; running median over +-0.25 x/D\n')
        for xv, zv in zip(xq_, np.interp(xq_, xx_, zz_)):
            fh.write(f'{xv:6.2f} {zv:8.4f}\n')
hdr = f'# x/D  z_x/D: first height (from the wall) where the conditional span-mean mu ratio RANS/193M = 1; source: {ZSRC}; running median +-{ZX_HW} D'
with open(os.path.join(DDIR, 'near_wall_boundary_zx.txt'), 'w') as fh:
    fh.write(hdr + '\n')
    for xv, zv in zip(*ZXTAB):
        fh.write(f'{xv:6.2f} {zv:8.4f}\n')
_xx, _zz = ZLINE
with open(os.path.join(DDIR, 'near_wall_boundary_zx_native.txt'), 'w') as fh:
    fh.write(hdr + ' (native x resolution)\n')
    for xv, zv in zip(_xx, _zz):
        if 3 <= xv <= 20 and np.isfinite(zv):
            fh.write(f'{xv:8.4f} {zv:8.4f}\n')
if MASK == '3d':
    with open(os.path.join(DDIR, 'near_wall_boundary_zx_pointfile.txt'), 'w') as fh:
        fh.write('# x/D  z_x/D from the point-file conditional span means (native x resolution)\n')
        for xv, zv in zip(XFP, ZXN):
            if 3 <= xv <= 20 and np.isfinite(zv):
                fh.write(f'{xv:8.4f} {zv:8.4f}\n')
