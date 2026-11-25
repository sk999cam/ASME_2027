#!/usr/bin/env python3
"""F19 fig_budget_lines: turbulence budget and length scales along wall-normal lines (Results 4.2).

Shows, against z/D (0 to 3), on the centreline (y/D = 0) and at mid-pitch (y/D = 3) at x/D = 4, 10, 20
(columns: the three stations; rows: centreline and mid-pitch for P_k and eps, then the same for the lengths;
panels lettered (a)-(l) row by row):
  (a)-(f)   production P_k (line) and dissipation eps (markers) of the turbulence kinetic energy, both per unit
            volume, normalised by rho_inf U_inf^3 / D (rho_inf = 0.900 kg/m3, U_inf = 86.2812 m/s, D = 7.75 mm);
  (g)-(l)   L_eff = k^1/2 |b| / |S*| (line; dns_fields.py of analysis/meff_c50b_c100b) and the dissipation length
            C_mu k^1.5 / eps (markers), C_mu = 0.09, both over D. Drawn at every point (no |S*| threshold), on a
            logarithmic axis because L_eff grows without bound where |S*| tends to zero towards the free stream.

Rungs (colour and style fixed by the common brief): 193M 3DNS (reference, stat3 batch b010, 10 flow-throughs), 50M 3DNS (b006),
20M LES C20bW (b006), RANS 50M (p0dns).  The 100M 3DNS and RANS 5M are not shown.
Complete: every rung is drawn from its own data; nothing is interpolated across rungs.

Definitions
  3DNS / LES  Favre k = R_ii/2 and P_k = -rho R_ij dU_i/dx_j (resolved stress, Favre mean velocity), eps =
              <Phi> - <tau_ij> d<u_i>/dx_j (dns/python/tke_budget.py, the solver's derivative operator on a
              sub-box around each line, margin 8 points, so the stencils equal the full-block ones).
              The stat3 sums are the cumulative mean over b006 (6 flow-throughs; 10 for the 193M, b010).
  20M LES     tau_ij carries the WALE stress, so eps is the TOTAL dissipation of the resolved k: resolved molecular
              part plus the sub-grid part (stat3 column 47 with the total tau, as the budget closes). The sub-grid share
              is estimated as eps_sgs = <Phi_sgs> (stat3 column 76) - 2 <mu_t> S*:S* and quoted in the report only.
              k is the resolved k (WALE has no sub-grid k); L_eff adds the mean WALE stress -2 <mu_t> S*, as dns_fields.py.
  RANS 50M    The RANS ran on the C50b mesh, so its vertex fields (rho, u, v, w, mu_t, k; analysis/meff_ransp0_c50b,
              rans_<n>.npy) go through the same sub-box procedure as the 3DNS (solver derivative operator).
              P_k = mu_t S^2 with S^2 = 2 S*_ij S*_ij; no production limiter and no gamma-Re_theta factor (the plain
              SST terms). omega comes from the STAR cell export (cells_rans_p0dns_full.csv, the same case), interpolated
              to each line point (log omega, inverse-distance weights of the 8 nearest cell centres); eps = beta* rho k
              omega, beta* = 0.09 (so C_mu k^1.5/eps = k^1/2/omega). L_eff of the RANS is the length behind its
              Boussinesq stress: mu_t = rho k^1/2 L_eff with cos theta = 1, so L_eff = mu_t / (rho k^1/2).
              The vertex k and mu_t and the cell omega are two discretisations of one field (agreement of
              L_eff with C_mu k^1.5/eps within 5 % in the film shows it).
Points on the mid-pitch line are taken at y/D = +-(3 - 0.01) and averaged (the pitch edge is a block face).
Check of the extraction: the 193M values at y/D = 0 differ from the cross-plane budget of dns/python/tke_budget_pf.py
(pictures_job14/stat3_b006/tke_budget/budget_planes.npz, x/D = 4 and 10) by at most 9 % in P_k and 2 % in eps and k;
that check was made with the b006 extraction of the 193M and was not repeated at b010.

Cache: /mnt/Elements/DNS/analysis/asme2027_figs/fig_budget_lines/ (npz per rung; the 193M file is lines_193M_b010.npz, the b006 file is kept); delete a file to recompute.
Run: python3 fig_budget_lines.py   (no arguments; the first run reads stat3 sub-boxes, minutes on the 193M).
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, '/home/sk999/Documents/4dns_filmcooled/dns/python')

D = 0.00775
U_INF = 86.2812
RHO_INF = 0.900                # free-stream density of the cases at z/D = 3 (0.8998-0.9011 kg/m3 in all rungs)
C_MU = 0.09
CACHE = '/mnt/Elements/DNS/analysis/asme2027_figs/fig_budget_lines'
OUT = '/home/sk999/Documents/ASME_2027/figures/fig_budget_lines.pdf'
XS = [4.0, 10.0, 20.0]
YL = {'cl': [0.0], 'mp': [3.0 - 0.01, -(3.0 - 0.01)]}
ZD = np.unique(np.concatenate([np.linspace(0.0, 0.3, 61), np.linspace(0.3, 3.2, 146)]))
RUNS = {'193M': ('/mnt/Elements/DNS/film_cooling/C193b_BR15', False),
        '50M': ('/mnt/Tank/DNS/film_cooling/C50b_BR15', False),
        '20M': ('/mnt/Tank/DNS/film_cooling/C20bW_BR15', True)}
BATCH_193M = 10                # reference batch of the 193M 3DNS (brief of 2026-11-25); the other rungs keep b006
BATCH = {'193M': BATCH_193M, '50M': 6, '20M': 6}
MARGIN = 8
COLS = ['rho', 'rhou', 'rhov', 'rhow', 'u', 'v', 'w', 'rhouu', 'rhovv', 'rhoww', 'rhouv', 'rhouw', 'rhovw',
        'phi', 'txx', 'tyy', 'tzz', 'txy', 'txz', 'tyz']
COLS_LES = ['mut', 'phi_sgs']


def line_points():
    """(n, 3) points in metres, and the line label (x, 'cl'|'mp', side) and z/D of each."""
    pts, lab, zz = [], [], []
    for xd in XS:
        for ln, ys in YL.items():
            for ys_ in ys:
                for z in ZD:
                    pts.append((xd * D, ys_ * D, z * D))
                    lab.append((xd, ln))
                    zz.append(z)
    return np.array(pts), lab, np.array(zz)


# ---------------------------------------------------------------- 3DNS / LES
def locate(run, label, pts, blocks):
    """Corners and weights of the points on the grid of run, cached (k-d tree over the whole grid once)."""
    path = os.path.join(CACHE, f'locate_{label}.npz')
    if os.path.exists(path):
        z = np.load(path)
        return z['cor'], z['w'], z['found']
    import plane_locate as pl
    # only the blocks whose corner bounding box (0.5 D margin) holds a line point: a k-d tree over all 193M
    # points needs more than the memory budget; the cells of the lines lie inside these blocks
    cand = []
    for b in blocks:
        gm = np.memmap(os.path.join(run, f'grid_{b["n"]}'), '<f8', 'r', shape=(b['nk'], b['nj'], b['ni'], 3))
        c8 = np.array([gm[k, j, i] for k in (0, -1) for j in (0, -1) for i in (0, -1)])
        lo, hi = c8.min(0) - 0.5 * D, c8.max(0) + 0.5 * D
        if np.any(np.all((pts >= lo) & (pts <= hi), axis=1)):
            cand.append(b)
        del gm
    print(f'{label}: candidate blocks {[b["n"] for b in cand]}', flush=True)
    gs = pl.GridSampler(run, cand)
    cs, w, found = gs.locate(pts)
    off_full = np.cumsum([0] + [b['nk'] * b['nj'] * b['ni'] for b in blocks])
    pos = {b['n']: q for q, b in enumerate(blocks)}
    bs = np.searchsorted(gs.off, cs, side='right') - 1
    cor = off_full[[pos[cand[i]['n']] for i in bs.ravel()]].reshape(cs.shape) + (cs - gs.off[bs])
    np.savez(path, cor=cor, w=w, found=found)
    del gs
    return cor, w, found


def fields_box(S, mt, X, les):
    """k, P, eps, L_eff and helpers on a sub-box. S: dict of column sums (nk, nj, ni)."""
    from fd_ops_f import metrics, grad
    m = {c: v / mt for c, v in S.items()}
    M = metrics(X)
    r = m['rho']
    U = np.stack([m['rhou'], m['rhov'], m['rhow']], -1) / r[..., None]
    ur = np.stack([m['u'], m['v'], m['w']], -1)
    rR = np.empty(r.shape + (3, 3))
    for (a, b), n in {(0, 0): 'rhouu', (1, 1): 'rhovv', (2, 2): 'rhoww', (0, 1): 'rhouv', (0, 2): 'rhouw',
                      (1, 2): 'rhovw'}.items():
        rR[..., a, b] = rR[..., b, a] = m[n] - r * U[..., a] * U[..., b]
    k = 0.5 * np.trace(rR, axis1=-2, axis2=-1) / r
    tau = np.empty(r.shape + (3, 3))
    for (a, b), n in {(0, 0): 'txx', (1, 1): 'tyy', (2, 2): 'tzz', (0, 1): 'txy', (0, 2): 'txz', (1, 2): 'tyz'}.items():
        tau[..., a, b] = tau[..., b, a] = m[n]
    gU = grad(U, M)                                      # [..., i, j] = dU_i/dx_j
    gur = grad(ur, M)
    P = -np.einsum('...ij,...ij', rR, gU)
    eps = m['phi'] - np.einsum('...ij,...ij', tau, gur)
    Sm = 0.5 * (gU + np.swapaxes(gU, -1, -2))
    tr = np.trace(Sm, axis1=-2, axis2=-1)
    for i in range(3):
        Sm[..., i, i] -= tr / 3
    a = rR.copy()
    for i in range(3):
        a[..., i, i] -= 2.0 / 3.0 * r * k
    out = {'rho': r, 'k': k, 'P': P, 'eps': eps}
    S2 = np.einsum('...ij,...ij', Sm, Sm)
    if les:
        mut = m['mut']
        a -= 2.0 * mut[..., None, None] * Sm
        out['eps_sgs'] = m['phi_sgs'] - 2.0 * mut * S2
        out['mut'] = mut
    an = np.sqrt(np.einsum('...ij,...ij', a, a))
    with np.errstate(divide='ignore', invalid='ignore'):
        out['L'] = np.sqrt(k) * (an / (2 * r * k)) / np.sqrt(S2)
        out['S'] = np.sqrt(S2)
    return out


RANS_VERT = '/mnt/Tank/DNS/analysis/meff_ransp0_c50b'     # rans_<n>.npy on the C50b grid (rans_on_grid.py)


def fields_rans(R, X):
    """RANS 50M on a sub-box of the C50b grid. R (nk, nj, ni, 8): rho, T, u, v, w, mu_t, k, mu.
    P = mu_t S^2 with S^2 = 2 S*_ij S*_ij (traceless strain, the solver's derivative operator)."""
    from fd_ops_f import metrics, grad
    M = metrics(X)
    U = np.ascontiguousarray(R[..., 2:5])
    G = grad(U, M)
    Sm = 0.5 * (G + np.swapaxes(G, -1, -2))
    tr = np.trace(Sm, axis1=-2, axis2=-1)
    for i in range(3):
        Sm[..., i, i] -= tr / 3
    S2 = np.einsum('...ij,...ij', Sm, Sm)
    rho, mut, k = R[..., 0], R[..., 5], R[..., 6]
    with np.errstate(divide='ignore', invalid='ignore'):
        L = mut / (rho * np.sqrt(k))                       # L_eff with cos = 1 (mu_eff = rho k^1/2 L_eff)
    return {'rho': rho, 'k': k, 'mut': mut, 'P': mut * 2.0 * S2, 'S': np.sqrt(S2), 'L': L}


def extract_dns(label, rans=False):
    """Line values of one rung (rans=True: the RANS 50M vertex fields on the C50b grid, label '50M')."""
    bt = BATCH[label]
    path = os.path.join(CACHE, 'lines_RANS_grid.npz' if rans else
                        (f'lines_{label}_b{bt:03d}.npz' if label == '193M' else f'lines_{label}.npz'))
    if os.path.exists(path):
        return dict(np.load(path))
    from tke_io_f import load_case, stat3_meta, stat3_ncol
    from tke_io import IDX
    run, les = RUNS[label]
    t0 = time.time()
    os.makedirs(CACHE, exist_ok=True)
    _, bd = load_case(run)
    blocks = [bd[n] for n in sorted(bd)]
    pts, lab, zz = line_points()
    cor, w, found = locate(run, label, pts, blocks)
    print(f'{label}: located {int(found.sum())}/{len(pts)} [{time.time() - t0:.0f} s]', flush=True)
    shp = [(b['nk'], b['nj'], b['ni']) for b in blocks]
    off = np.cumsum([0] + [int(np.prod(s)) for s in shp])
    names = (['rho', 'k', 'mut', 'P', 'L', 'S'] if rans else
             ['rho', 'k', 'P', 'eps', 'L', 'S'] + (['eps_sgs', 'mut'] if les else []))
    val = {n: np.full(len(pts), np.nan) for n in names}
    cols = COLS + (COLS_LES if les else [])
    # group points by (line, block of corner 0)
    cc = cor[:, 0]
    blk = np.searchsorted(off, cc, side='right') - 1
    lines = {}
    for p in np.nonzero(found)[0]:
        lines.setdefault((lab[p], blk[p], p // len(ZD)), []).append(p)
    for (ln, bi, _), idx in sorted(lines.items(), key=lambda t: (t[0][1], t[0][2])):
        idx = np.array(idx)
        b = blocks[bi]
        nk, nj, ni = shp[bi]
        loc = cor[idx] - off[bi]
        kk, r_ = np.divmod(loc, nj * ni)
        jj, ii = np.divmod(r_, ni)
        if np.any((cor[idx] < off[bi]) | (cor[idx] >= off[bi + 1])):
            raise RuntimeError('cell straddles blocks')
        lo = [max(0, a.min() - MARGIN) for a in (kk, jj, ii)]
        hi = [min(n, a.max() + MARGIN + 1) for a, n in zip((kk, jj, ii), (nk, nj, ni))]
        if hi[0] - lo[0] < 10 or hi[1] - lo[1] < 10 or hi[2] - lo[2] < 10:
            raise RuntimeError('sub-box too small')
        sl = np.s_[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
        gm = np.memmap(os.path.join(run, f'grid_{b["n"]}'), '<f8', 'r', shape=(nk, nj, ni, 3))
        X = np.array(gm[sl])
        if rans:
            R = np.load(os.path.join(RANS_VERT, f'rans_{b["n"]}.npy'), mmap_mode='r')
            f = fields_rans(np.array(R[sl]), X)
        else:
            ncol = stat3_ncol(run, b, f'_b{bt:03d}')
            mm = np.memmap(os.path.join(run, f'stat3_{b["n"]}_b{bt:03d}'), '<f8', 'r', shape=(nk, nj, ni, ncol))
            sub = np.array(mm[sl])
            S = {c: np.ascontiguousarray(sub[..., IDX[c]]) for c in cols}
            mt = stat3_meta(run, b['n'], f'_b{bt:03d}')['mean_time']
            del sub, mm
            f = fields_box(S, mt, X, les)
        del gm
        sh = f['k'].shape
        for n in names:
            flat = f[n].ravel()
            loc_ = ((kk[:, :] - lo[0]) * sh[1] + (jj - lo[1])) * sh[2] + (ii - lo[2])
            v8 = flat[loc_]
            val[n][idx] = (v8 * w[idx]).sum(axis=1)
        print(f'{label}: block {b["n"]} line {ln} {len(idx)} points sub-box {tuple(h - l for h, l in zip(hi, lo))} '
              f'[{time.time() - t0:.0f} s]', flush=True)
    keys = [f'{xd:g}_{ln}' for xd in XS for ln in YL]
    out = {'zD': ZD, 'found': found.reshape(len(XS), -1)}
    # average the two mid-pitch sides, then arrange (station, line, z)
    for n in names:
        a = val[n].reshape(len(XS), -1)               # per station: cl (nz), mp+ (nz), mp- (nz)
        nz = len(ZD)
        cl = a[:, :nz]
        mp = np.nanmean(np.stack([a[:, nz:2 * nz], a[:, 2 * nz:3 * nz]]), axis=0)
        out[n] = np.stack([cl, mp], axis=1)            # (station, line, z)
    np.savez(path, **out)
    return out


# ---------------------------------------------------------------- RANS
RANS_CSV = '/mnt/Tank/STARCCM/film_cooling/rans_50M_Lm25b_BR15_FST_ZPG_p0dns/cells_rans_p0dns_full.csv'


def rans_cells():
    """Cells of the STAR export in boxes around the six lines, cached."""
    path = os.path.join(CACHE, 'rans_cells.npz')
    if os.path.exists(path):
        return dict(np.load(path))
    import pyarrow.csv as pacsv
    os.makedirs(CACHE, exist_ok=True)
    keys = {'Velocity[i]': 'u', 'Velocity[j]': 'v', 'Velocity[k]': 'w', 'Temperature': 'T',
            'Turbulent Kinetic Energy': 'k', 'Specific Dissipation Rate': 'om', 'Turbulent Viscosity': 'mut',
            'Intermittency': 'g', 'Density': 'rho', 'Pressure': 'p', 'Volume': 'V', 'X (m)': 'x', 'Y (m)': 'y',
            'Z (m)': 'z'}
    rd = pacsv.open_csv(RANS_CSV, read_options=pacsv.ReadOptions(block_size=1 << 27))
    names = rd.schema.names
    col = {}
    for kq, v in keys.items():
        h = [n for n in names if n.startswith(kq)]
        assert len(h) == 1, (kq, names)
        col[v] = h[0]
    out = {v: [] for v in keys.values()}
    t0, n = time.time(), 0
    for batch in rd:
        x = batch.column(col['x']).to_numpy() / D
        y = batch.column(col['y']).to_numpy() / D
        z = batch.column(col['z']).to_numpy() / D
        m = (z < 3.8)
        mx = np.zeros(len(x), bool)
        for xs in XS:
            mx |= np.abs(x - xs) < 0.7
        m &= mx & ((np.abs(y) < 0.7) | (np.abs(y) > 2.3))
        if m.any():
            for v, c in col.items():
                out[v].append(batch.column(c).to_numpy()[m])
        n += len(x)
        if time.time() - t0 > 60 * getattr(rans_cells, 'k', 1):
            rans_cells.k = getattr(rans_cells, 'k', 1) + 1
            print(f'RANS {n / 1e6:.0f}M rows [{time.time() - t0:.0f} s]', flush=True)
    out = {v: np.concatenate(a) for v, a in out.items()}
    np.savez(path, **out)
    print(f'RANS cells kept {len(out["x"])} of {n} [{time.time() - t0:.0f} s]')
    return out


def rans_lines(nn=8):
    """RANS 50M on the six lines. k, mu_t, rho, S from the vertex fields on the C50b grid (same operator as the
    3DNS); omega and gamma from the STAR cell export, interpolated to each line point (log omega, inverse-distance
    weights of the nn nearest cell centres). eps = beta* rho k omega."""
    path = os.path.join(CACHE, 'lines_RANS.npz')
    if os.path.exists(path):
        return dict(np.load(path))
    from scipy.spatial import cKDTree
    out = extract_dns('50M', rans=True)
    c = rans_cells()
    xyz = np.stack([c['x'], c['y'], c['z']], 1) / D
    tree = cKDTree(xyz)
    nz = len(ZD)
    om = np.full((len(XS), 2, nz), np.nan)
    g = np.full((len(XS), 2, nz), np.nan)
    for a, xd in enumerate(XS):
        for b, (ln, ys) in enumerate(YL.items()):
            oo, gg = [], []
            for ysd in ys:
                P = np.stack([np.full(nz, xd), np.full(nz, ysd), ZD], 1)
                d, ix = tree.query(P, k=nn)
                w = 1.0 / (d + 1e-6) ** 2
                w /= w.sum(axis=1, keepdims=True)
                oo.append(np.exp((w * np.log(c['om'][ix])).sum(axis=1)))
                gg.append((w * c['g'][ix]).sum(axis=1))
            om[a, b] = np.mean(oo, axis=0)
            g[a, b] = np.mean(gg, axis=0)
    out['om'], out['g'] = om, g
    out['eps'] = C_MU * out['rho'] * out['k'] * om
    np.savez(path, **out)
    return out


# ---------------------------------------------------------------- figure
RUNG = [('193M', '193M 3DNS', '#000000', '-'), ('50M', '50M 3DNS', '#D55E00', ':'),
        ('20M', '20M LES', '#E69F00', '-.'), ('RANS', 'RANS 50M', '#009E73', '-')]


def norm_terms(d, key):
    """P_k, eps (both over rho_inf U^3/D), L_eff / D and C_mu k^1.5/eps / D of one rung, arrays (station, line, z).
    No threshold: the lengths are drawn at every point (non-finite values, if any, are not plotted)."""
    unit = RHO_INF * U_INF ** 3 / D
    rho = d['rho'] if 'rho' in d else RHO_INF
    P = d['P'] / unit
    eps = d['eps'] / unit
    with np.errstate(divide='ignore', invalid='ignore'):
        L = d['L'] / D
        Leps = C_MU * d['k'] ** 1.5 / (d['eps'] / rho) / D
    L = np.where(np.isfinite(L) & (L > 0), L, np.nan)
    Leps = np.where(np.isfinite(Leps) & (Leps > 0), Leps, np.nan)
    return P, eps, L, Leps


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from cycler import cycler
    sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
    from functions import set_mpl_params
    data = {k: extract_dns(k) for k in ('193M', '50M', '20M')}
    data['RANS'] = rans_lines()
    set_mpl_params(0.4, grid_on=True)
    plt.rcParams['axes.prop_cycle'] = cycler(color=['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442',
                                                    '#0072B2', '#D55E00', '#CC79A7'])
    fig, axs = plt.subplots(4, 3, figsize=(7.2, 8.2), sharex=True)
    tag = [f'({chr(97 + i)})' for i in range(12)]      # (a)-(l), row by row
    rowinfo = [(0, 'P', r'$y/D=0$'), (1, 'P', r'$y/D=3$'), (0, 'L', r'$y/D=0$'), (1, 'L', r'$y/D=3$')]
    for r, (line, kind, ltxt) in enumerate(rowinfo):
        for c, xd in enumerate(XS):
            ax = axs[r, c]
            for key, lab, col, ls in RUNG:
                d = data[key]
                P, eps, L, Le = norm_terms(d, key)
                z = d['zD']
                m = z <= 3.0
                a, b = (1e3 * P, 1e3 * eps) if kind == 'P' else (L, Le)
                ax.plot(z[m], a[c, line][m], color=col, ls=ls)
                ax.plot(z[m], b[c, line][m], color=col, ls=ls, marker='o', markevery=15, markerfacecolor='none')
            if r == 0:
                ax.text(0.5, 1.03, rf'$x/D={xd:g}$', transform=ax.transAxes, ha='center', va='bottom')
            ax.text(0.96, 0.94, f'{tag[3 * r + c]} {ltxt}', transform=ax.transAxes, ha='right', va='top')
            ax.set_xlim(0, 3)
            if kind == 'L':
                ax.set_yscale('log')
                ax.set_ylim(1e-2, 3)     # all points are drawn; max drawn L is about 1.1 D, the wall values below 1e-2 lie outside the view
            if c == 0:
                ax.set_ylabel(r'$10^3\,(P_k,\ \varepsilon)\,D/(\rho_\infty U_\infty^3)$' if kind == 'P'
                              else r'$(L_\mathrm{eff},\ C_\mu k^{3/2}/\varepsilon)\,/D$')
            if r == 3:
                ax.set_xlabel(r'$z/D$')
    h = [plt.Line2D([], [], color=col, ls=ls, label=lab) for _, lab, col, ls in RUNG]
    h += [plt.Line2D([], [], color='k', ls='-', label=r'line: $P_k$, $L_\mathrm{eff}$'),
          plt.Line2D([], [], color='k', ls='none', marker='o', markerfacecolor='none',
                     label=r'markers: $\varepsilon$, $C_\mu k^{3/2}/\varepsilon$')]
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    fig.legend(handles=h, loc='lower center', ncol=3, bbox_to_anchor=(0.5, 0.0))
    fig.savefig(OUT)
    print('saved', OUT)


if __name__ == '__main__':
    main()
