"""20M LES counterparts of the RANS 50M / 193M numbers quoted in the manuscript (table only, no figure).

Usage: python3 les_counterparts_table.py      -> figures/code/data/les_counterparts_table.txt

For every quantity: the manuscript's RANS and 193M values, the same-route reproduction of them (route check), the LES value, the sources.
Routes are those of the scripts named in each section (fig_mean_features.py, n7_series files, fig_wall_state.py, fig_meff_split_film.py,
fig_prt_data tables, fig_cross_planes.py). Nothing is approximated: where the LES data do not exist the text says what is needed.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
import matplotlib  # noqa: E402
matplotlib.use('Agg')
import fig_coolant_budget as CB  # noqa: E402
import _meff_span_split as SS2  # noqa: E402

MF = CB.MF
AS = '/mnt/Elements/DNS/analysis/asme2027_figs/'
TAB = os.path.join(HERE, 'data', 'les_counterparts_table.txt')
L = []
P = L.append


def rng(a, m):
    return f'{np.nanmin(a[m]):.2f} to {np.nanmax(a[m]):.2f}'


# ------------------------------------------------------------------------------------------------ mean-flow features (a, c, d)
out = MF.run()
lvol, lpl = CB.les_analysis()
XS = MF.XS
q = {'RANS 50M': out['rans']['q'], '193M': out['dns']['q'], '20M LES': lvol['q']}
m420 = (XS >= 4 - 1e-9) & (XS <= 20 + 1e-9)
ix = lambda x: int(np.argmin(np.abs(XS - x)))      # noqa: E731

P('# 20M LES counterparts of the manuscript quantities (figures/code/les_counterparts_table.py). RANS = RANS 50M, 193M = C193b b010, LES = C20bW (b006).')
P('# Each block: manuscript values, same-route reproduction (route check), LES value, sources. "LES data not available" = needs data from the data session.')
P('')
P('=' * 110)
P('a. CRVP core height (manuscript: the RANS cores lie ~0.1 D below those of the reference, not shown)')
P('   Route: fig_mean_features.py (analyse): omega_x, Gaussian SIG = (2, 3, 3) cells, vorticity-weighted centroid inside 0.5 of the peak, |y| < 1.5, z >= 0.05;')
P('   mean of the two sides; n4 grid (0.05 D). LES: same functions on mean_features/les20_mean_grid.npz (ResampleToImage of the b006 block means).')
P(f'   {"x/D":>5s} {"zc 193M":>9s} {"zc RANS":>9s} {"zc LES":>9s} {"RANS-193M":>10s} {"LES-193M":>9s} {"RANS/193M":>10s} {"LES/193M":>9s}')
for x in (4, 6, 10, 15, 20):
    a, r, l_ = (q[k]['zc/D'][ix(x)] for k in ('193M', 'RANS 50M', '20M LES'))
    P(f'   {x:5d} {a:9.4f} {r:9.4f} {l_:9.4f} {r - a:+10.3f} {l_ - a:+9.3f} {r / a:10.3f} {l_ / a:9.3f}')
d_r = q['RANS 50M']['zc/D'] - q['193M']['zc/D']
d_l = q['20M LES']['zc/D'] - q['193M']['zc/D']
mm = (XS >= 6 - 1e-9) & (XS <= 20 + 1e-9)
P(f'   mean of (core height - 193M) over x/D 6-20: RANS {np.nanmean(d_r[mm]):+.3f} D (manuscript: ~0.1 D below, reproduced), LES {np.nanmean(d_l[mm]):+.3f} D')
P('   Sources: RANS/193M data/mean_features_table.txt (zc/D rows); LES mean_features/les20_mean_grid.npz.')
P('')
P('=' * 110)
P('c. Centreline upwash at the core height and inward velocity at z/D 0.1 under the pair (manuscript: RANS 1.14-1.66 and 1.75-1.93 times the reference over x/D 4-20)')
P('   Route: fig_mean_features.cross_to: W/U on the centreline y = 0 interpolated at the mean core height zc; V/U inward = mean over the sides of the maximum inward V/U at z/D 0.1, |y| < 1.5.')
for key, nm in (('W/U centre at zc', 'upwash W/U at zc'), ('V/U inward z=0.1', 'inward V/U at z/D 0.1')):
    rr, ll = q['RANS 50M'][key] / q['193M'][key], q['20M LES'][key] / q['193M'][key]
    st = [ix(x) for x in (4, 6, 10, 15, 20)]
    P(f'   {nm}: ratio to 193M over the stations x/D 4, 6, 10, 15, 20: RANS {np.min(rr[st]):.2f} to {np.max(rr[st]):.2f}; LES {np.min(ll[st]):.2f} to {np.max(ll[st]):.2f}'
      f' (all steps of 0.25 D over x/D 4-20: RANS {rng(rr, m420)}, LES {rng(ll, m420)})')
    P(f'     {"x/D":>5s} {"193M":>9s} {"RANS":>9s} {"LES":>9s} {"RANS/193M":>10s} {"LES/193M":>9s}')
    for x in (4, 6, 10, 15, 20):
        a, r, l_ = (q[k][key][ix(x)] for k in ('193M', 'RANS 50M', '20M LES'))
        P(f'     {x:5d} {a:9.4f} {r:9.4f} {l_:9.4f} {r / a:10.3f} {l_ / a:9.3f}')
P('   Sources: data/mean_features_table.txt (RANS, 193M); LES: mean_features/les20_mean_grid.npz.')
P('')
P('=' * 110)
P('d. Plume half-width and coolant-centroid height (manuscript: at x/D 20 the RANS half-width is 0.75 of the reference, the centroid height agrees within 12 %)')
P('   Route: fig_mean_features.cross_to: half-width = lateral half width at 0.5 of the centreline value of the column integral of theta (clipped at 0), mean of the two sides;')
P('   centroid = theta-weighted mean height of the plane.')
P(f'   {"x/D":>5s} | {"half width: 193M":>16s} {"RANS":>8s} {"LES":>8s} {"RANS/193":>9s} {"LES/193":>8s} | {"centroid: 193M":>14s} {"RANS":>8s} {"LES":>8s} {"RANS/193":>9s} {"LES/193":>8s}')
for x in (4, 6, 10, 15, 20):
    h = [q[k]['half width'][ix(x)] for k in ('193M', 'RANS 50M', '20M LES')]
    c_ = [q[k]['z centroid'][ix(x)] for k in ('193M', 'RANS 50M', '20M LES')]
    P(f'   {x:5d} | {h[0]:16.4f} {h[1]:8.4f} {h[2]:8.4f} {h[1] / h[0]:9.3f} {h[2] / h[0]:8.3f} | {c_[0]:14.4f} {c_[1]:8.4f} {c_[2]:8.4f} {c_[1] / c_[0]:9.3f} {c_[2] / c_[0]:8.3f}')
P('   Sources: data/mean_features_table.txt (RANS, 193M); LES: mean_features/les20_mean_grid.npz (theta).')

# ------------------------------------------------------------------------------------------------ b. fluxes through z/D planes
P('')
P('=' * 110)
P('b. Vertical coolant flux through the planes z/D 0.25, 0.5, 1.0 (manuscript: RANS turbulent flux 0.67-0.99 of the reference through z/D 1.0, towards the wall through z/D 0.25;')
P('   ordering of the two parts held through the three planes)')
P('   Route: wall_balance/n7_series_c193b_b010_rans50.npz (keys h<h>_193_Jm/Jt, h<h>_rans_Jm/Jt; n7_summary.txt), g s^-1 m^-1, one pitch, x/D 2-20 every 0.25 D. LES: wall_balance/')
P('   n7_series_c193b_b010_les20.npz (keys h<h>_les20_Jm, _Jt, _Jsgs and _se, six-batch jackknife of the LES; turbulent = resolved Jt + mean SGS Jsgs, the SGS part has no SE;')
P('   the SE quoted for the LES turbulent flux is that of the resolved part). Same stations as the 193M and RANS series, no interpolation.')
s = np.load(AS + 'wall_balance/n7_series_c193b_b010_rans50.npz')
ls_ = np.load(AS + 'wall_balance/n7_series_c193b_b010_les20.npz')
xs = s['x_D']
assert np.allclose(xs, ls_['x_D'])
mx = (xs >= 4 - 1e-9) & (xs <= 20 + 1e-9)
for h in ('0.25', '0.5', '1.0'):
    j193m, j193t = s[f'h{h}_193_Jm'], s[f'h{h}_193_Jt']
    jrm, jrt = s[f'h{h}_rans_Jm'], s[f'h{h}_rans_Jt']
    jlm = ls_[f'h{h}_les20_Jm']
    jlt = ls_[f'h{h}_les20_Jt'] + ls_[f'h{h}_les20_Jsgs']
    jlm_i, jlt_i = jlm, jlt
    se193m, se193t = s[f'h{h}_193_Jm_se'], s[f'h{h}_193_Jt_se']
    selm, selt = ls_[f'h{h}_les20_Jm_se'], ls_[f'h{h}_les20_Jt_se']
    first = lambda a, b: (f'{xs[mx][np.where(a[mx] > b[mx])[0][0]]:.2f}' if (a[mx] > b[mx]).any() else 'not in 4-20')      # noqa: E731
    P(f'   z/D = {h}: ratio to 193M over x/D 4-20 (steps 0.25): turbulent flux RANS {rng(jrt / j193t, mx)}, LES {rng(jlt_i / j193t, mx)}; '
      f'mean-flow flux RANS {rng(jrm / j193m, mx)}, LES {rng(jlm_i / j193m, mx)}')
    P(f'      first x/D (>= 4) where turbulent > mean flux: 193M {first(j193t, j193m)}, RANS {first(jrt, jrm)}, LES {first(jlt_i, jlm_i)}; '
      f'RANS turbulent flux negative (towards the wall) at {int((jrt[mx] < 0).sum())} of {int(mx.sum())} stations, LES {int((jlt_i[mx] < 0).sum())}')
    P(f'      LES turbulent flux negative at {int((jlt[mx] < 0).sum())} of {int(mx.sum())} stations; LES mean flux negative at {int((jlm[mx] < 0).sum())} (first at x/D '
      f'{xs[mx][np.where(jlm[mx] < 0)[0][0]]:.2f})' if (jlm[mx] < 0).any() else f'      LES turbulent flux negative at {int((jlt[mx] < 0).sum())} of {int(mx.sum())} stations; LES mean flux never negative in x/D 4-20')
    P(f'      {"x/D":>5s} | {"193M Jm (SE)":>16s} {"193M Jt (SE)":>16s} | {"RANS Jm":>8s} {"RANS Jt":>8s} | {"LES Jm (SE)":>16s} {"LES Jt (SE)":>16s}')
    for x in (5, 10, 20):
        i = int(np.argmin(np.abs(xs - x)))
        P(f'      {x:5d} | {j193m[i]:8.2f} ({se193m[i]:5.2f}) {j193t[i]:8.2f} ({se193t[i]:5.2f}) | {jrm[i]:8.2f} {jrt[i]:8.2f} | '
          f'{jlm_i[i]:8.2f} ({selm[i]:5.2f}) {jlt_i[i]:8.2f} ({selt[i]:5.2f})')
P('   Surface z = z_x: see coolant_budget_les_table.txt (six-batch LES, with SE).')

# ------------------------------------------------------------------------------------------------ e. gamma
P('')
P('=' * 110)
P('e. Plane-mean intermittency gamma at z/D 0.1, x/D -3 to 20, |y/D| < 3 (manuscript: reference 0.73, RANS 0.17; reference 0.726, RANS 0.175 as quoted by the coordinator)')
import _gamma_extract_v as GV  # noqa: E402
import fig_wall_state as WS  # noqa: E402
gx = np.arange(-3.0, 20.0 + 1e-9, 0.04)
gy = np.arange(-3.0, 3.0 + 1e-9, 0.04)
r = WS.rans_slab()
sel = (r['z'] > 0.094) & (r['z'] < 0.102)
r = {k: v[sel] for k, v in r.items()}
g_r = np.nanmean(WS.regrid(r['x'], r['y'], r['g'], gx, gy, 0.15))
res = {}
for tag in ('c193b', 'c20bw'):
    G = GV.load_z010(tag)['gamma']
    i1 = np.clip(np.round((gx + 3.0) / 0.02).astype(int), 0, GV.NXR - 1)
    i2 = np.clip(np.round((gy + 3.0) / 0.02).astype(int), 0, GV.NY - 1)
    res[tag] = (float(np.nanmean(G[np.ix_(i1, i2)])), float(G[GV.XG <= 20 + 1e-9].mean()))
P('   Route (fig_wall_state.py main): gamma maps on the grid x/D -3..20 and |y/D| <= 3 in steps of 0.04 D, plain mean. RANS: Intermittency (SST gamma-Re_theta model variable) of the prism')
P('   layer z/D 0.098 (cells_rans_p0dns_full.csv through rans50_z010_cells.npz), regridded with fig_wall_state.regrid (mask distance 0.15 D). 193M and LES: high-pass spanwise-stress')
P('   detector, Otsu threshold pooled over the 193M band |y/D| 2-3, x/D 0-20 (one threshold for both), z010 plane sampler (_gamma_extract_v.load_z010).')
P(f'   Reproduction: RANS {g_r:.4f} (manuscript 0.17; the value 0.175 is NOT reproduced: variants of the regrid give 0.167-0.180, mask distance 0.05 D gives 0.180),')
P(f'   193M {res["c193b"][0]:.4f} on the 0.04 D grid (manuscript 0.73) and {res["c193b"][1]:.4f} on the native 0.02 D detector grid (x <= 20; the quoted 0.726)')
P(f'   LES: {res["c20bw"][0]:.4f} on the same grid (native detector grid: {res["c20bw"][1]:.4f}); LES - 193M = {res["c20bw"][0] - res["c193b"][0]:+.4f}')
P('   Sources: gamma_v/c193b_z010_v.npz, gamma_v/c20bw_z010_v.npz, gamma_v/thresholds.json; RANS: bl_extract/rans50_z010_cells.npz.')

# ------------------------------------------------------------------------------------------------ f. eddy-viscosity factor split
P('')
P('=' * 110)
P('f. Eddy-viscosity factor split, ratios to the 193M of k^1/2, L (L_eff) and 1/cos(theta) below and above z_x (fig_meff_split_film panels (b), (c) give them for the RANS)')
import _meff_split_nearfar as SP  # noqa: E402
zt = np.loadtxt(os.path.join(HERE, 'data', 'near_wall_boundary_zx.txt'))
ZXT = (zt[:, 0], zt[:, 1])
FMIN = 0.2
CY = dict(mu='cy_mu_mean', rho='cy_rho_mean', k='cy_sqrtk_mean', L='cy_Leff_mean', c='cy_cos_mean')
k0 = dict(SS2.KEYS)
SS2.KEYS = CY
try:
    GR, GL = SS2.load(case='RANS50'), SS2.load(case='C20bW')
finally:
    SS2.KEYS = k0


def span_ratios(G, lay, xc, width):
    """Ratios (run / 193M) of the layer means with the weights of _meff_span_split.compute_curve (area x valid fraction, f >= FMIN) and the mask of its split."""
    XC, ZC = np.meshgrid(G['xc'], G['zc'])
    zx = np.interp(XC, *ZXT)
    ok = np.isfinite(G['f']) & (G['f'] >= FMIN) & (ZC > 0) & (ZC <= 2 + 1e-9)
    for k in SS2.KEYS:
        ok &= np.isfinite(G['R'][k]) & np.isfinite(G['D'][k])
    ok &= (G['R']['mu'] > 0) & (G['D']['mu'] > 0) & (G['R']['rho'] > 0) & (G['D']['rho'] > 0) & (G['R']['k'] > 0) & (G['D']['k'] > 0)
    ok &= (G['R']['L'] > 0) & (G['D']['L'] > 0) & (G['D']['c'] > 0)
    slab = ((ZC < zx) if lay == 'near' else (ZC >= zx)) & (np.abs(XC - xc) < width / 2) & (ZC > 0) & (ZC <= 2 + 1e-9) & np.isfinite(G['f'])
    m = slab & ok
    w = (G['area'] * np.where(np.isfinite(G['f']), G['f'], 0.0))[m]
    mean = lambda a: (w * a[m]).sum() / w.sum()      # noqa: E731
    return dict(k=mean(G['R']['k']) / mean(G['D']['k']), L=mean(G['R']['L']) / mean(G['D']['L']), icos=mean(G['R']['c']) / mean(G['D']['c']),
                mu=mean(G['R']['mu']) / mean(G['D']['mu']))


P('   Routes: RANS (paper, panels b, c): point file, 3D strain mask, mesh-volume weights (_meff_split_nearfar.compute_curve, vmean: k^1/2, L, 1/cos = 1/<cos_193>, RANS cos = 1).')
P('   RANS and LES by the span route: meff_span files, Delta-y weighted conditional means (cy), rebinned to 0.1 D, weights dx dz f (f >= 0.2), ratio of the layer means; the 1/cos column is <cos_run>/<cos_193> (= 1/<cos_193> for the RANS)')
P('   (RANS cos = 1; LES cos is its own: cos = -a:S*/(|a||S*|) of the resolved stress plus mean WALE stress as in the README of meff_span). The point-file route has no LES counterpart.')
P('   Layer means over x/D 4-20 (one slab 16 D wide) and slabs of 1 D at x/D 5, 10, 19.5; layers: below z_x (wall to z_x), above z_x (to z/D 2 where f >= 0.2).')
P(f'   {"layer":>6s} {"x/D":>8s} | {"RANS point file: k1/2":>21s} {"L":>6s} {"1/cos":>6s} {"mu":>6s} | {"RANS span: k1/2":>16s} {"L":>6s} {"1/cos":>6s} {"mu":>6s} | {"LES span: k1/2":>15s} {"L":>6s} {"1/cos":>6s} {"mu":>6s}')
for lay in ('near', 'above'):
    for lab, xc, wd in (('4-20', 12.0, 16.0), ('5', 5.0, 1.0), ('10', 10.0, 1.0), ('19.5', 19.5, 1.0)):
        pf = SP.compute_curve(ZXT, [xc], wd, layer=lay)[0]
        vm = pf['vmean']
        a = (vm['k'][0] / vm['k'][1], vm['L'][0] / vm['L'][1], 1.0 / vm['c'][1], float(np.exp(pf['tot'])))
        rs, ls = span_ratios(GR, lay, xc, wd), span_ratios(GL, lay, xc, wd)
        P(f'   {lay:>6s} {lab:>8s} | {a[0]:21.3f} {a[1]:6.3f} {a[2]:6.3f} {a[3]:6.3f} | {rs["k"]:16.3f} {rs["L"]:6.3f} {rs["icos"]:6.3f} {rs["mu"]:6.3f} | '
          f'{ls["k"]:15.3f} {ls["L"]:6.3f} {ls["icos"]:6.3f} {ls["mu"]:6.3f}')
cur = SP.compute_curve(ZXT, np.arange(4.5, 19.5 + 1e-9, 0.25), 1.0, layer='near')
cura = SP.compute_curve(ZXT, np.arange(4.5, 19.5 + 1e-9, 0.25), 1.0, layer='above')
rk = lambda c_, f_: np.array([f_(r_['vmean']) for r_ in c_])      # noqa: E731
P(f'   Panel (b)/(c) ranges of the RANS over the slabs x/D 4.5-19.5 (point file): below z_x k^1/2 {rk(cur, lambda v: v["k"][0] / v["k"][1]).min():.2f}-{rk(cur, lambda v: v["k"][0] / v["k"][1]).max():.2f}, '
  f'L {rk(cur, lambda v: v["L"][0] / v["L"][1]).min():.2f}-{rk(cur, lambda v: v["L"][0] / v["L"][1]).max():.2f}, 1/cos {rk(cur, lambda v: 1 / v["c"][1]).min():.2f}-{rk(cur, lambda v: 1 / v["c"][1]).max():.2f}; '
  f'above z_x k^1/2 {rk(cura, lambda v: v["k"][0] / v["k"][1]).min():.2f}-{rk(cura, lambda v: v["k"][0] / v["k"][1]).max():.2f}, L {rk(cura, lambda v: v["L"][0] / v["L"][1]).min():.2f}-{rk(cura, lambda v: v["L"][0] / v["L"][1]).max():.2f}, '
  f'1/cos {rk(cura, lambda v: 1 / v["c"][1]).min():.2f}-{rk(cura, lambda v: 1 / v["c"][1]).max():.2f}')
lr = {lay: {k_: [] for k_ in ('k', 'L', 'icos')} for lay in ('near', 'above')}
for lay in lr:
    for xc in np.arange(4.5, 19.5 + 1e-9, 0.25):
        s_ = span_ratios(GL, lay, xc, 1.0)
        for k_ in lr[lay]:
            lr[lay][k_].append(s_[k_])
P('   LES (span route) ranges over the same slabs: ' + '; '.join(f'{lay} ' + ', '.join(f'{k_} {min(v):.2f}-{max(v):.2f}' for k_, v in lr[lay].items()) for lay in lr))
P('   Sources: RANS point file: fig_meff_split_film.py (pts.npz via _meff_split_nearfar.py); span: meff_span/span_RANS50.npz, span_C20bW.npz, span_C193b.npz; z_x: data/near_wall_boundary_zx.txt.')

# ------------------------------------------------------------------------------------------------ g, h
P('')
P('=' * 110)
P('g. Cross-plane turbulence level and length scale (manuscript, near L706: on the x/D 4, 10, 20 planes the reference spreads near-wall turbulence from the jet to mid-pitch, the RANS has')
P('   little turbulence outside the jet, the short-length band along the wall is thicker in the RANS; no number is quoted)')
import meff_planes_common as mc  # noqa: E402
DAT = {'RANS 50M': mc.load('rans'), '193M': mc.load('ref'), '20M LES': np.load(AS + 'meff_planes/les20_b006.npz', allow_pickle=True)}


def fld(qn, rung, p):
    dd = DAT[rung]
    v = dd[f'{p}_{qn}'].astype(float)
    if qn == 'k':
        return v / mc.U_INF ** 2
    return np.where((dd[f'{p}_k'] > 0) & (v > 0), v / mc.D, np.nan)


P('   Route (fig_cross_planes.py print): planewise median and 99th percentile of k/U_inf^2 and L_eff/D (RANS: mu_t/(rho k^1/2)) over the plane points (z/D 0-3).')
P('   LES: meff_planes/les20_b006.npz (C20bW stat3 b006, same keys and plane points as rans50_p0dns.npz; k and L_eff as for mu_eff with the mean WALE stress added), same field() as fig_cross_planes.')
P(f'   {"x/D":>5s} {"qty":>4s} | {"RANS med":>10s} {"RANS p99":>10s} | {"193M med":>10s} {"193M p99":>10s} | {"LES med":>10s} {"LES p99":>10s} | {"RANS/193M med":>14s} {"p99":>6s} | {"LES/193M med":>13s} {"p99":>6s}')
for x in (4, 10, 20):
    for qn in ('k', 'L'):
        a, b, c_ = fld(qn, 'RANS 50M', f'x{x}'), fld(qn, '193M', f'x{x}'), fld(qn, '20M LES', f'x{x}')
        P(f'   {x:5d} {qn:>4s} | {np.nanmedian(a):10.4g} {np.nanpercentile(a, 99):10.4g} | {np.nanmedian(b):10.4g} {np.nanpercentile(b, 99):10.4g} | {np.nanmedian(c_):10.4g} {np.nanpercentile(c_, 99):10.4g} | '
          f'{np.nanmedian(a) / np.nanmedian(b):14.3f} {np.nanpercentile(a, 99) / np.nanpercentile(b, 99):6.3f} | {np.nanmedian(c_) / np.nanmedian(b):13.3f} {np.nanpercentile(c_, 99) / np.nanpercentile(b, 99):6.3f}')
P('')
P('=' * 110)
P('h. Spanwise diffusivities (manuscript, near L746: alpha_y/alpha_z of the reference is 3.2 on the film flanks at x/D 10; nu_y/nu_z 0.7-1.1; the lateral heat diffusivity of the RANS at the film')
P('   edges is 0.5-0.6 times the reference)')
P('   Route: fig_prt_data.py (regional least squares of -rho v\'\'T\'\' dT/dy etc. on the cross-planes x/D 10 and 20, flanks |y|/D >= 1.5, z/D >= 0.3, b010); table prt_classes/prt_data_tables.txt.')
P('   Manuscript values reproduced from that table: alpha_y/alpha_z flanks x/D 10 = 3.16 (3.2), nu_y/nu_z flanks 1.07 (x/D 10) and 0.71 (x/D 20) (0.7-1.1);')
P('   RANS/193M lateral heat diffusivity at the edges = exp(ln(P_R/P_D) - gradient term) = exp(0.041 - 0.658) = 0.54 (x/D 10), exp(0.191 - 0.650) = 0.63 (x/D 20) (0.5-0.6).')
import prt_corr_common as cc  # noqa: E402
import fig_prt_classes as fc  # noqa: E402
Rr = {'193M': dict(np.load(f'{mc.CACHE}/ref193_b010.npz')), '20M LES': dict(np.load(f'{mc.CACHE}/les20_b006.npz', allow_pickle=True))}
RANSd = dict(np.load(f'{mc.CACHE}/rans50_p0dns.npz'))
IX = {k: i for i, k in enumerate(fc.KEYS)}
NUR = (IX['lnm'], IX['ldm'], IX['cnm'], IX['cdm'])         # nu_y/nu_z, as fig_prt_data
ALR = (IX['lna'], IX['lda'], IX['cna'], IX['cda'])         # alpha_y/alpha_z
ratio = lambda S, spec: (S[spec[0]].sum() / S[spec[1]].sum()) / (S[spec[2]].sum() / S[spec[3]].sum())      # noqa: E731
rg = np.random.default_rng(5)


def boot(S, spec, nb=1000):
    have = np.nonzero(S[IX['cdm']] > 0)[0]
    Sh = S[:, have]
    w = rg.multinomial(len(have), np.ones(len(have)) / len(have), size=nb).astype(float)
    v = ((w @ Sh[spec[0]]) / (w @ Sh[spec[1]])) / ((w @ Sh[spec[2]]) / (w @ Sh[spec[3]]))
    v = v[np.isfinite(v)]
    return np.percentile(v, [2.5, 97.5])


def load_corr_tag(p, tag, shape):
    """as fig_prt_classes.load_corr, with the file name of the LES (corr_<plane>_c20bw_b006.npz)."""
    z = dict(np.load(os.path.join(cc.CORR_DIR, tag.format(p=p))))
    out = {}
    for qq in cc.FIELDS:
        a = np.full(int(np.prod(shape)), np.nan, np.float32)
        a[z['idx']] = z[qq]
        out[qq] = a.astype(float)
    return out


FILES = {'193M': 'corr_{p}_b010.npz', '20M LES': 'corr_{p}_c20bw_b006.npz'}
EDGE = lambda d: (d['z'] >= 0.3) & (np.abs(d['x']) >= 1.0)      # noqa: E731
CL = [('all, z/D >= 0.3', lambda d: d['z'] >= 0.3), ('centre |y|/D < 1', lambda d: (d['z'] >= 0.3) & (np.abs(d['x']) < 1)),
      ('flanks |y|/D >= 1.5', lambda d: (d['z'] >= 0.3) & (np.abs(d['x']) >= 1.5)), ('wall z/D < 0.3', lambda d: d['z'] < 0.3)]
P('   Route (fig_prt_data.py, functions of fig_prt_classes.py: prep, sums): per-plane sums of -rho u\'\'v\'\' dU/dy, (dU/dy)^2, -rho u\'\'w\'\' dU/dz, (dU/dz)^2 and the heat analogues over the')
P('   points of a class, ratio of the sums; 95 % interval from the block bootstrap over 0.5 D y bins (1000 draws, seed 5; the same code for both runs). LES: corr_x10/x20_c20bw_b006.npz')
P('   (SGS parts added to the fluxes and shear stresses) with meff_planes/les20_b006.npz where the 193M uses ref193_b010.npz; 193M: corr_x10/x20_b010.npz.')
P('   Route differences forced by the LES files: none in the formulas. The LES is one batch (b006, 6 FT) against the 193M b010 (10 FT); intervals are bootstrap over y bins for both (no')
P('   jackknife over batches exists for either); the LES planes are the same 600 x 300 trilinear samples of the (coarser) 20M mesh, so neighbouring plane points are more strongly')
P('   correlated than for the 193M and the bootstrap interval of the LES is optimistic by an unquantified factor (mesh spacing against plane spacing not checked here).')
res = {}
for p in ('x10', 'x20'):
    for nm in ('193M', '20M LES'):
        Rn = Rr[nm]
        Pd = fc.prep(Rn, load_corr_tag(p, FILES[nm], Rn[f'{p}_A'].shape), p)
        res[(p, nm)] = Pd
P(f'   {"plane":>5s} {"class":>20s} | {"nu_y/nu_z 193M":>15s} {"[95 %]":>15s} {"LES":>7s} {"[95 %]":>15s} | {"alpha_y/alpha_z 193M":>20s} {"[95 %]":>15s} {"LES":>7s} {"[95 %]":>15s}')
for p in ('x10', 'x20'):
    for lab, f in CL:
        row = []
        for nm in ('193M', '20M LES'):
            S, n = fc.sums(res[(p, nm)], f(res[(p, nm)]))
            row.append((ratio(S, NUR), boot(S, NUR), ratio(S, ALR), boot(S, ALR)))
        a, b = row
        P(f'   {p:>5s} {lab:>20s} | {a[0]:15.2f} [{a[1][0]:5.2f},{a[1][1]:6.2f}] {b[0]:7.2f} [{b[1][0]:5.2f},{b[1][1]:6.2f}] | {a[2]:20.2f} [{a[3][0]:5.2f},{a[3][1]:6.2f}] {b[2]:7.2f} [{b[3][0]:5.2f},{b[3][1]:6.2f}]')
P('   Reproduction of the manuscript values from this code (193M): alpha_y/alpha_z flanks x/D 10 and nu_y/nu_z flanks x/D 10, 20 are the 193M entries of the flanks rows above (manuscript 3.2; 0.7-1.1;')
P('   table prt_classes/prt_data_tables.txt: 3.16, 1.07, 0.71).')
P('')
P('   Lateral heat diffusivity at the film edges (|y|/D >= 1, z/D >= 0.3): alpha_y = sum(-rho v\'\'T\'\' dT/dy) / sum (dT/dy)^2 over the edge points; RANS: sum(mu_t/0.9 (dT/dy)^2) / sum (dT/dy)^2')
P('   with dT/dy the finite difference of the RANS plane temperature (fig_prt_data part b; ratio RANS/193M = exp(ln(P_R/P_D) - gradient term)); LES: same sums as the 193M from its correlations.')
P(f'   {"plane":>5s} | {"RANS/193M (manuscript 0.5-0.6)":>31s} | {"LES/193M":>9s} | {"alpha_y 193M":>13s} {"RANS":>11s} {"LES":>11s} (kg m^-1 s^-1 per K-gradient sums; mu_t/0.9 for RANS)')
for p in ('x10', 'x20'):
    d = res[(p, '193M')]
    m = EDGE(d) & d['ok']
    S, n = fc.sums(d, EDGE(d))
    aD = S[IX['lna']].sum() / S[IX['lda']].sum()
    T = RANSd[f'{p}_T']
    dy = RANSd[f'{p}_A'][0, 1] - RANSd[f'{p}_A'][0, 0]
    gR = np.gradient(T, dy, axis=1).ravel()
    mut = RANSd[f'{p}_meff'].ravel()
    aR = np.sum((mut[m] / 0.9) * gR[m] ** 2) / np.sum(gR[m] ** 2)
    dl = res[(p, '20M LES')]
    Sl, nl = fc.sums(dl, EDGE(dl))
    aL = Sl[IX['lna']].sum() / Sl[IX['lda']].sum()
    P(f'   {p:>5s} | {aR / aD:31.3f} | {aL / aD:9.3f} | {aD:13.3e} {aR:11.3e} {aL:11.3e}   (points: 193M {n}, LES {nl})')
P('   Reproduction: RANS/193M 0.54 (x/D 10) and 0.63 (x/D 20) from the table of the manuscript route; compare the RANS/193M column.')
P('')
P('=' * 110)
P('Items needing new data from the data session: none of a-h is left open. Limits: b LES series SE excludes the SGS part; h LES is one batch (b006) and the bootstrap interval is optimistic (see above);')
P('f has the span route only for the LES (the point-file route has no LES counterpart); e: the RANS 0.175 quoted by the coordinator is not reproduced (0.167 by the fig_wall_state route).')
open(TAB, 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
print('saved', TAB)
