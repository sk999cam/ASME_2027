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
    P(f'   {nm}: ratio to 193M over x/D 4-20: RANS {rng(rr, m420)} (steps of 0.25 D); LES {rng(ll, m420)}')
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
P('   Route: wall_balance/n7_series_c193b_b010_rans50.npz (keys h<h>_193_Jm/Jt, h<h>_rans_Jm/Jt; n7_summary.txt), g s^-1 m^-1, one pitch. LES: wall_balance/les_c50b/')
P('   balance_c20bw_b1.npz and _b2.npz (keys h<h>_Jm, h<h>_Jt, h<h>_Jsgs, kg s^-1 m^-1 x 1000; mean of the two batches; turbulent = resolved + mean SGS).')
P('   LIMIT: only two of the LES batches are available as plane series (batch-to-batch scatter of Jm 10-40 %, e.g. h 0.5, x/D 5: 11.7 and 10.2); the six-batch LES')
P('   series (SE, like the 193M ten-batch jackknife) needs the remaining batches from the data session. Values below are two-batch means without SE.')
s = np.load(AS + 'wall_balance/n7_series_c193b_b010_rans50.npz')
b1 = np.load(AS + 'wall_balance/les_c50b/balance_c20bw_b1.npz')
b2 = np.load(AS + 'wall_balance/les_c50b/balance_c20bw_b2.npz')
xs, xl = s['x_D'], b1['xs']
mx = (xs >= 4 - 1e-9) & (xs <= 20 + 1e-9)
for h in ('0.25', '0.5', '1.0'):
    j193m, j193t = s[f'h{h}_193_Jm'], s[f'h{h}_193_Jt']
    jrm, jrt = s[f'h{h}_rans_Jm'], s[f'h{h}_rans_Jt']
    jlm = 1e3 * 0.5 * (b1[f'h{h}_Jm'] + b2[f'h{h}_Jm'])
    jlt = 1e3 * 0.5 * (b1[f'h{h}_Jt'] + b2[f'h{h}_Jt'] + b1[f'h{h}_Jsgs'] + b2[f'h{h}_Jsgs'])
    jlm_i, jlt_i = np.interp(xs, xl, jlm), np.interp(xs, xl, jlt)
    first = lambda a, b: (f'{xs[mx][np.where(a[mx] > b[mx])[0][0]]:.2f}' if (a[mx] > b[mx]).any() else 'not in 4-20')      # noqa: E731
    P(f'   z/D = {h}: ratio to 193M over x/D 4-20 (steps 0.25): turbulent flux RANS {rng(jrt / j193t, mx)}, LES {rng(jlt_i / j193t, mx)}; '
      f'mean-flow flux RANS {rng(jrm / j193m, mx)}, LES {rng(jlm_i / j193m, mx)}')
    P(f'      first x/D (>= 4) where turbulent > mean flux: 193M {first(j193t, j193m)}, RANS {first(jrt, jrm)}, LES {first(jlt_i, jlm_i)}; '
      f'RANS turbulent flux negative (towards the wall) at {int((jrt[mx] < 0).sum())} of {int(mx.sum())} stations, LES {int((jlt_i[mx] < 0).sum())}')
    P(f'      {"x/D":>5s} {"193M Jm":>8s} {"Jt":>7s} {"RANS Jm":>8s} {"Jt":>7s} {"LES Jm":>8s} {"Jt":>7s}')
    for x in (5, 10, 20):
        i = int(np.argmin(np.abs(xs - x)))
        P(f'      {x:5d} {j193m[i]:8.2f} {j193t[i]:7.2f} {jrm[i]:8.2f} {jrt[i]:7.2f} {jlm_i[i]:8.2f} {jlt_i[i]:7.2f}')
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
P(f'   Reproduction: RANS {g_r:.4f} (manuscript 0.17), 193M {res["c193b"][0]:.4f} (manuscript 0.73 / 0.726); also on the native detector grid (x <= 20, 0.02 D): 193M {res["c193b"][1]:.4f}')
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
    return dict(k=mean(G['R']['k']) / mean(G['D']['k']), L=mean(G['R']['L']) / mean(G['D']['L']), icos=mean(G['D']['c']) / mean(G['R']['c']),
                mu=mean(G['R']['mu']) / mean(G['D']['mu']))


P('   Routes: RANS (paper, panels b, c): point file, 3D strain mask, mesh-volume weights (_meff_split_nearfar.compute_curve, vmean: k^1/2, L, 1/cos = 1/<cos_193>, RANS cos = 1).')
P('   RANS and LES by the span route: meff_span files, Delta-y weighted conditional means (cy), rebinned to 0.1 D, weights dx dz f (f >= 0.2), ratio of the layer means; 1/cos = <cos_193>/<cos_run>')
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
DAT = {'RANS 50M': mc.load('rans'), '193M': mc.load('ref')}


def fld(qn, rung, p):
    dd = DAT[rung]
    v = dd[f'{p}_{qn}'].astype(float)
    if qn == 'k':
        return v / mc.U_INF ** 2
    return np.where((dd[f'{p}_k'] > 0) & (v > 0), v / mc.D, np.nan)


P('   Route (fig_cross_planes.py print): planewise median and 99th percentile of k/U_inf^2 and L_eff/D (RANS: mu_t/(rho k^1/2)) over the plane points (z/D 0-3).')
P(f'   {"x/D":>5s} {"quantity":>9s} | {"RANS median":>12s} {"RANS p99":>10s} | {"193M median":>12s} {"193M p99":>10s} | {"RANS/193M median":>17s}')
for x in (4, 10, 20):
    for qn in ('k', 'L'):
        a, b = fld(qn, 'RANS 50M', f'x{x}'), fld(qn, '193M', f'x{x}')
        P(f'   {x:5d} {qn:>9s} | {np.nanmedian(a):12.4g} {np.nanpercentile(a, 99):10.4g} | {np.nanmedian(b):12.4g} {np.nanpercentile(b, 99):10.4g} | {np.nanmedian(a) / np.nanmedian(b):17.3f}')
P('   LES data not available: needs k and L_eff = k^1/2 |b|/|S*| (stat3 correlations of C20bW b006, with the WALE stress added to the stress as for mu_eff) on the x/D 4, 10, 20 planes in the')
P('   format of meff_planes/ref193_b010.npz (keys x4_k, x4_L, x4_S, x4_T, x4_A, x4_B ...) from the data session. les20_mean_grid.npz has no k or Reynolds stress.')
P('')
P('=' * 110)
P('h. Spanwise diffusivities (manuscript, near L746: alpha_y/alpha_z of the reference is 3.2 on the film flanks at x/D 10; nu_y/nu_z 0.7-1.1; the lateral heat diffusivity of the RANS at the film')
P('   edges is 0.5-0.6 times the reference)')
P('   Route: fig_prt_data.py (regional least squares of -rho v\'\'T\'\' dT/dy etc. on the cross-planes x/D 10 and 20, flanks |y|/D >= 1.5, z/D >= 0.3, b010); table prt_classes/prt_data_tables.txt.')
P('   Manuscript values reproduced from that table: alpha_y/alpha_z flanks x/D 10 = 3.16 (3.2), nu_y/nu_z flanks 1.07 (x/D 10) and 0.71 (x/D 20) (0.7-1.1);')
P('   RANS/193M lateral heat diffusivity at the edges = exp(ln(P_R/P_D) - gradient term) = exp(0.041 - 0.658) = 0.54 (x/D 10), exp(0.191 - 0.650) = 0.63 (x/D 20) (0.5-0.6).')
P('   LES data not available: needs the plane correlations of fig_prt_data / prt_corr_common (sum of -rho u\'\'v\'\' dU/dy, -rho u\'\'w\'\' dU/dz, -rho v\'\'T\'\' dT/dy, -rho w\'\'T\'\' dT/dz, with the')
P('   mean gradients) for C20bW stat3 b006 on the planes x/D 10 and 20, with the SGS stress and heat flux added; no corr_*_C20bW file exists in prt_classes/.')
P('')
P('=' * 110)
P('Items needing new data: b (the full six-batch LES plane series with SE; two batches exist), g (LES k and L_eff on the cross-planes x/D 4, 10, 20), h (LES plane correlations x/D 10, 20).')
P('Items computed for the LES: a, c, d (mean grid), e (gamma z010 planes), f (meff_span, span route only: the point-file route has no LES counterpart), b at reduced statistics.')
open(TAB, 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
print('saved', TAB)
