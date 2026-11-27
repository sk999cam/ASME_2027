#!/usr/bin/env python3
"""CANDIDATE fig_budget_maps: maps of turbulence production P_k and dissipation eps (Results 4.2), replacement
candidate for the 12-panel line figure fig_budget_lines (which is not touched).

Usage
  python3 fig_budget_maps.py plan            cost of the 193M read (no data read), prints sub-box sizes
  python3 fig_budget_maps.py rans            RANS 50M extraction (Tank, read-only), cached
  python3 fig_budget_maps.py ref             193M 3DNS (stat3 b010) extraction, cached per block, resumable
  python3 fig_budget_maps.py check           quadrature and extraction checks against caches (prints)
  python3 fig_budget_maps.py [A|B|AB]        figure(s) (VARIANT; default AB), medians printed
Outputs: figures/fig_budget_maps.png (variant A, log10(P_k/eps), 6 panels)
         figures/fig_budget_maps_B.png (variant B, P_k and eps on one log scale, 12 panels).

Planes (x/D -2 to 25, z/D 0 to 2, uniform 0.05 D x 0.01 D grid; pitch 6 D, hole centre at y/D = 0)
  centreline  y/D = 0
  mid-pitch   mean of y/D = +2.99 and -2.99 (the pitch edge y/D = +-3 is a block face; as fig_budget_lines)
  span-avg    trapezoid rule over y/D at the levels  0, +-1, +-2, +-2.99  (7 planes, SEL), divided by the span
              covered (5.98 D). A sampled average, not the full-volume average: the 193M stat3 sub-volume of the
              whole pitch is 48 GB, the 7 strips 17.6 GB. The same levels are used for both rungs (the RANS cache
              holds 13 levels, 0.5 D apart, of which SEL is used).
Quantities, normalised by rho_inf U_inf^3 / D (rho_inf = 0.900 kg/m3, U_inf = 86.2812 m/s, D = 7.75 mm)
  193M 3DNS  Favre k, P_k = -rho R_ij dU_i/dx_j (resolved stress, Favre mean velocity), eps = <Phi> - <tau_ij>
             d<u_i>/dx_j, stat3 batch b010 (10 flow-throughs), solver derivative operator on sub-boxes with a
             margin of 8 points (fig_budget_lines.fields_box, unchanged). Buffer points (mean tau = 0): NaN.
  RANS 50M   P_k = mu_t * 2 S*_ij S*_ij (S* the traceless mean strain rate; plain SST, no limiter, no gamma-Re_theta
             factor), eps = beta* rho k omega, beta* = C_mu = 0.09 (fig_budget_lines.fields_rans; omega from the
             STAR cell export cells_rans_p0dns_full.csv, log-omega inverse-distance weights of the 8 nearest cells).
Ratio P_k/eps for span-average = <P_k>/<eps> (ratio of the averages).
Negligible cells (|P_k| and eps both below THRESH, absolute, all panels) are drawn flat light grey in variant A
(the ratio is meaningless there); the data stay in the cache and in variant B.
Cache: /mnt/Elements/DNS/analysis/asme2027_figs/budget_maps/
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, '/home/sk999/Documents/4dns_filmcooled/dns/python')
import fig_budget_lines as fbl      # noqa: E402  (definitions and field operators reused unchanged)

D = fbl.D
U_INF, RHO_INF, C_MU = fbl.U_INF, fbl.RHO_INF, fbl.C_MU
UNIT = RHO_INF * U_INF ** 3 / D
CACHE = '/mnt/Elements/DNS/analysis/asme2027_figs/budget_maps'
EPS_CACHE = '/mnt/Tank/DNS/analysis/eps_prt_c193b/b010'        # full-volume eps (W/m3, float32) of the 193M, b010
RUN193 = fbl.RUNS['193M'][0]
RUN50 = fbl.RUNS['50M'][0]
BATCH = 10
MARGIN = 8
OUTDIR = '/home/sk999/Documents/ASME_2027/figures'
X0, X1, DX = -2.0, 25.0, 0.05
Z0, Z1, DZ = 0.0, 2.0, 0.01
XG = np.round(np.arange(X0, X1 + 1e-9, DX), 6)
ZG = np.round(np.arange(Z0, Z1 + 1e-9, DZ), 6)
YL = np.array([-2.99, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 2.99])
ICL = 6
# levels used for the span average and read from the 193M stat3: 13 levels would cost 30.4 GB, 9 levels 22.7 GB,
# these 7 (y/D = -2.99, -2, -1, 0, 1, 2, 2.99) cost 17.6 GB. Trapezoid error of eps against the 13-level
# estimate (eps full-volume cache): median 1.6 %, 90th percentile 3.8 %.
SEL = [0, 2, 4, 6, 8, 10, 12]
# absolute threshold, same for every panel: a cell is negligible when |P_k| and eps are both below THRESH
# (units rho_inf U_inf^3 / D). 1e-5 (lower limit of figure B) greys almost nothing (193M eps >= 2e-5), so 3e-5 is used.
THRESH = 3e-5
NX, NZ, NL = len(XG), len(ZG), len(YL)
XPLOT = 20.0                       # domain ends at x/D 20 (buffer beyond): plotted and quoted range
XP = XG[XG <= XPLOT + 1e-9]
NP = len(XP)


def points():
    """(NL*NZ*NX, 3) points in metres, level-major, then z, then x."""
    Y, Z, X = np.meshgrid(YL, ZG, XG, indexing='ij')
    return np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1) * D


# ---------------------------------------------------------------- location (candidate blocks only)
def locate(run, label, pts, blocks):
    path = os.path.join(CACHE, f'locate_{label}.npz')
    if os.path.exists(path):
        z = np.load(path)
        return z['cor'], z['w'], z['found']
    import plane_locate as pl
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
    os.makedirs(CACHE, exist_ok=True)
    np.savez(path, cor=cor, w=w, found=found)
    del gs
    return cor, w, found


def merge_boxes(boxes):
    """Merge overlapping index boxes [(lo(3), hi(3))] into their bounding boxes (repeat until disjoint)."""
    boxes = [(np.array(l), np.array(h)) for l, h in boxes]
    done = False
    while not done:
        done = True
        out = []
        for lo, hi in boxes:
            for q, (l2, h2) in enumerate(out):
                if np.all(lo < h2) and np.all(l2 < hi):
                    out[q] = (np.minimum(lo, l2), np.maximum(hi, h2))
                    done = False
                    break
            else:
                out.append((lo, hi))
        boxes = out
    return boxes


def plan_boxes(cor, found, blocks, shp, off):
    """{block index: [(lo, hi, point indices)]}: merged sub-boxes (margin included) per block."""
    cc = cor[:, 0]
    blk = np.searchsorted(off, cc, side='right') - 1
    npl = NZ * NX
    plan = {}
    for bi in np.unique(blk[found]):
        idx = np.nonzero(found & (blk == bi))[0]
        nk, nj, ni = shp[bi]
        loc = cor[idx] - off[bi]
        if np.any((cor[idx] < off[bi]) | (cor[idx] >= off[bi + 1])):
            raise RuntimeError('cell straddles blocks')
        kk, r_ = np.divmod(loc, nj * ni)
        jj, ii = np.divmod(r_, ni)
        lev = idx // npl
        bx, plv = [], []
        for l in np.unique(lev):
            m = lev == l
            lo = [max(0, a[m].min() - MARGIN) for a in (kk, jj, ii)]
            hi = [min(n, a[m].max() + MARGIN + 1) for a, n in zip((kk, jj, ii), (nk, nj, ni))]
            bx.append((lo, hi))
        mb = merge_boxes(bx)
        res = []
        for lo, hi in mb:
            inb = np.all((kk >= lo[0] + 0) & (kk < hi[0]) & (jj >= lo[1]) & (jj < hi[1]) & (ii >= lo[2]) & (ii < hi[2]),
                         axis=1)
            res.append((lo, hi, idx[inb], kk[inb], jj[inb], ii[inb]))
        plan[bi] = res
    return plan


# ---------------------------------------------------------------- 193M
def run_ref(plan_only=False):
    from tke_io_f import load_case, stat3_ncol, stat3_meta
    from tke_io import IDX
    path = os.path.join(CACHE, f'ref193_b{BATCH:03d}_levels.npz')
    if os.path.exists(path) and not plan_only:
        return dict(np.load(path))
    t0 = time.time()
    os.makedirs(CACHE, exist_ok=True)
    _, bd = load_case(RUN193)
    blocks = [bd[n] for n in sorted(bd)]
    pts = points()
    cor, w, found = locate(RUN193, '193M', pts, blocks)
    print(f'193M: located {int(found.sum())}/{len(pts)} [{time.time() - t0:.0f} s]', flush=True)
    selmask = np.zeros(NL, bool)
    selmask[SEL] = True                      # read only the levels of SEL (cost limit)
    found = found & np.repeat(selmask, NZ * NX)
    shp = [(b['nk'], b['nj'], b['ni']) for b in blocks]
    off = np.cumsum([0] + [int(np.prod(s)) for s in shp])
    plan = plan_boxes(cor, found, blocks, shp, off)
    tot = 0
    for bi, res in plan.items():
        for lo, hi, *_ in res:
            tot += int(np.prod(np.array(hi) - np.array(lo)))
    print(f'193M plan: {sum(len(r) for r in plan.values())} sub-boxes, {tot / 1e6:.1f} M points, '
          f'{tot * 592 / 1e9:.1f} GB of stat3', flush=True)
    if plan_only:
        return None
    ncl = fbl.COLS
    pdir = os.path.join(CACHE, 'ref_parts')
    os.makedirs(pdir, exist_ok=True)
    val = {n: np.full(len(pts), np.nan) for n in ('P', 'eps')}
    for bi, res in sorted(plan.items()):
        b = blocks[bi]
        pf = os.path.join(pdir, f'blk_{b["n"]}.npz')
        if os.path.exists(pf):
            z = np.load(pf)
            for n in val:
                val[n][z['idx']] = z[n]
            print(f'block {b["n"]}: from part cache', flush=True)
            continue
        nk, nj, ni = shp[bi]
        mt = stat3_meta(RUN193, b['n'], f'_b{BATCH:03d}')['mean_time']
        ncol = stat3_ncol(RUN193, b, f'_b{BATCH:03d}')
        mm = np.memmap(os.path.join(RUN193, f'stat3_{b["n"]}_b{BATCH:03d}'), '<f8', 'r', shape=(nk, nj, ni, ncol))
        gm = np.memmap(os.path.join(RUN193, f'grid_{b["n"]}'), '<f8', 'r', shape=(nk, nj, ni, 3))
        bidx, bP, bE = [], [], []
        for lo, hi, idx, kk, jj, ii in res:
            sl = np.s_[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
            t1 = time.time()
            sub = np.array(mm[sl])
            S = {c: np.ascontiguousarray(sub[..., IDX[c]]) for c in ncl}
            buf = np.all(np.stack([sub[..., IDX[c]] for c in ('txx', 'tyy', 'tzz', 'txy', 'txz', 'tyz')], -1) == 0.0, -1)
            del sub
            X = np.array(gm[sl])
            f = fbl.fields_box(S, mt, X, False)
            del S, X
            sh = f['k'].shape
            loc_ = ((kk - lo[0]) * sh[1] + (jj - lo[1])) * sh[2] + (ii - lo[2])
            res_ = {}
            for n in ('P', 'eps'):
                fl = np.where(buf, np.nan, f[n]).ravel()
                res_[n] = (fl[loc_] * w[idx]).sum(axis=1)
            bidx.append(idx)
            bP.append(res_['P'])
            bE.append(res_['eps'])
            print(f'block {b["n"]}: sub-box {tuple(h - l for h, l in zip(hi, lo))} {len(idx)} points '
                  f'[{time.time() - t1:.0f} s, total {time.time() - t0:.0f} s]', flush=True)
            del f, buf
        idx = np.concatenate(bidx)
        np.savez(pf, idx=idx, P=np.concatenate(bP), eps=np.concatenate(bE))
        val['P'][idx] = np.concatenate(bP)
        val['eps'][idx] = np.concatenate(bE)
        del mm, gm
    out = {n: v.reshape(NL, NZ, NX) for n, v in val.items()}
    out['found'] = found.reshape(NL, NZ, NX)
    out['yl'], out['xg'], out['zg'], out['mean_time_batch'] = YL, XG, ZG, BATCH
    np.savez(path, **out)
    print(f'193M done [{time.time() - t0:.0f} s]', flush=True)
    return out


# ---------------------------------------------------------------- RANS 50M
def rans_cells_region():
    path = os.path.join(CACHE, 'rans_cells_region.npz')
    if os.path.exists(path):
        return dict(np.load(path))
    import pyarrow.csv as pacsv
    os.makedirs(CACHE, exist_ok=True)
    keys = {'Specific Dissipation Rate': 'om', 'X (m)': 'x', 'Y (m)': 'y', 'Z (m)': 'z'}
    rd = pacsv.open_csv(fbl.RANS_CSV, read_options=pacsv.ReadOptions(block_size=1 << 27))
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
        m = (z < 2.4) & (x > X0 - 0.4) & (x < X1 + 0.4) & (np.abs(y) < 3.4)
        if m.any():
            for v, c in col.items():
                out[v].append(batch.column(c).to_numpy()[m])
        n += len(x)
    out = {v: np.concatenate(a) for v, a in out.items()}
    np.savez(path, **out)
    print(f'RANS cells kept {len(out["x"])} of {n} [{time.time() - t0:.0f} s]', flush=True)
    return out


def run_rans():
    from tke_io_f import load_case
    from scipy.spatial import cKDTree
    path = os.path.join(CACHE, 'rans50_levels.npz')
    if os.path.exists(path):
        return dict(np.load(path))
    t0 = time.time()
    os.makedirs(CACHE, exist_ok=True)
    _, bd = load_case(RUN50)
    blocks = [bd[n] for n in sorted(bd)]
    pts = points()
    cor, w, found = locate(RUN50, '50M', pts, blocks)
    print(f'RANS: located {int(found.sum())}/{len(pts)} [{time.time() - t0:.0f} s]', flush=True)
    shp = [(b['nk'], b['nj'], b['ni']) for b in blocks]
    off = np.cumsum([0] + [int(np.prod(s)) for s in shp])
    plan = plan_boxes(cor, found, blocks, shp, off)
    val = {n: np.full(len(pts), np.nan) for n in ('rho', 'k', 'mut', 'P')}
    for bi, res in sorted(plan.items()):
        b = blocks[bi]
        nk, nj, ni = shp[bi]
        gm = np.memmap(os.path.join(RUN50, f'grid_{b["n"]}'), '<f8', 'r', shape=(nk, nj, ni, 3))
        R = np.load(os.path.join(fbl.RANS_VERT, f'rans_{b["n"]}.npy'), mmap_mode='r')
        for lo, hi, idx, kk, jj, ii in res:
            sl = np.s_[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
            f = fbl.fields_rans(np.array(R[sl]), np.array(gm[sl]))
            sh = f['k'].shape
            loc_ = ((kk - lo[0]) * sh[1] + (jj - lo[1])) * sh[2] + (ii - lo[2])
            for n in val:
                val[n][idx] = (f[n].ravel()[loc_] * w[idx]).sum(axis=1)
        print(f'RANS block {b["n"]} [{time.time() - t0:.0f} s]', flush=True)
    c = rans_cells_region()
    tree = cKDTree(np.stack([c['x'], c['y'], c['z']], 1) / D)
    d, ix = tree.query(pts / D, k=8, workers=-1)
    wt = 1.0 / (d + 1e-6) ** 2
    wt /= wt.sum(axis=1, keepdims=True)
    om = np.exp((wt * np.log(c['om'][ix])).sum(axis=1))
    out = {n: v.reshape(NL, NZ, NX) for n, v in val.items()}
    out['om'] = om.reshape(NL, NZ, NX)
    out['eps'] = C_MU * out['rho'] * out['k'] * out['om']
    out['found'] = found.reshape(NL, NZ, NX)
    out['yl'], out['xg'], out['zg'] = YL, XG, ZG
    np.savez(path, **out)
    print(f'RANS done [{time.time() - t0:.0f} s]', flush=True)
    return out


# ---------------------------------------------------------------- planes
def planes(d):
    """{'cl','mp','sa'} -> (P, eps) normalised, arrays (NZ, NX). Span average by the trapezoid rule over YL."""
    P, E = d['P'][..., :NP] / UNIT, d['eps'][..., :NP] / UNIT
    out = {'cl': (P[ICL], E[ICL]), 'mp': (0.5 * (P[0] + P[-1]), 0.5 * (E[0] + E[-1]))}
    span = YL[-1] - YL[0]
    y = YL[SEL]
    out['sa'] = (np.trapz(P[SEL], y, axis=0) / span, np.trapz(E[SEL], y, axis=0) / span)
    return out


def negligible(P, E):
    with np.errstate(invalid='ignore'):
        return (np.abs(P) < THRESH) & (np.abs(E) < THRESH)


def medians(pl):
    """Median P/eps in the film (x/D 3-20, 0.1 <= z/D <= 1.0, non-negligible cells), and in 3-10 and 10-20."""
    P, E = pl
    neg = negligible(P, E)
    X, Z = np.meshgrid(XP, ZG)
    base = (Z >= 0.1) & (Z <= 1.0) & ~neg & np.isfinite(P) & np.isfinite(E) & (E > 0)
    res = {}
    for nm, (a, b) in {'3-20': (3, 20), '3-10': (3, 10), '10-20': (10, 20)}.items():
        m = base & (X >= a) & (X <= b)
        res[nm] = float(np.median(P[m] / E[m])) if m.any() else np.nan
    res['frac_negl'] = float(np.mean(neg[(Z >= 0.1) & (Z <= 1.0) & (X >= 3)]))
    return res


# ---------------------------------------------------------------- checks
def check():
    ref, rs = run_ref(), run_rans()
    lines = dict(np.load('/mnt/Elements/DNS/analysis/asme2027_figs/fig_budget_lines/lines_193M_b010.npz'))
    print('found fraction 193M', ref['found'].mean(), 'RANS', rs['found'].mean())
    for si, xd in enumerate(fbl.XS):
        ix = int(round((xd - X0) / DX))
        # y=0 line of fig_budget_lines versus the plane at y = 0 (z/D 0.05..2)
        zl = lines['zD']
        for nm in ('P', 'eps'):
            a = np.interp(ZG, zl, lines[nm][si, 0])
            b = ref[nm][ICL][:, ix]
            m = (ZG > 0.05) & (ZG < 2) & (np.abs(a) > 0.02 * np.nanmax(np.abs(a)))
            print(f'x/D {xd:g}: 193M {nm} plane vs lines cache, max rel diff {np.nanmax(np.abs(a[m] - b[m]) / np.abs(a[m])):.3f}')
    # quadrature error of the span average, eps, against the full-volume cache
    from tke_io_f import load_case
    _, bd = load_case(RUN193)
    print('eps full-volume cache check: see check_eps_quadrature()')


# ---------------------------------------------------------------- figure
OKABE = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
ROWS = [('cl', r'$y/D=0$'), ('mp', r'$y/D=\pm3$'), ('sa', r'span average, $|y/D|\leq3$')]
HOLE = (-1.0, 2.36)       # hole exit on the centre plane, z = 0 (fig_config_regions EXIT polygon)


def setup():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from cycler import cycler
    sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
    from functions import set_mpl_params
    set_mpl_params(0.4, grid_on=False)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE)
    return plt


def mark_hole(ax, row):
    if row == 'cl':
        ax.plot(HOLE, [0, 0], color='k', lw=2.2, solid_capstyle='butt', clip_on=False)
    ax.set_xlim(X0, XPLOT)
    ax.set_ylim(Z0, Z1)


def figure_A(plt, pls):
    from matplotlib.colors import ListedColormap, Normalize
    fig, axs = plt.subplots(6, 1, figsize=(7.2, 6.4), sharex=True)
    rungs = [('rans', 'RANS 50M'), ('ref', '193M 3DNS')]
    k = 0
    for r, (row, rl) in enumerate(ROWS):
        for key, name in rungs:
            ax = axs[k]
            P, E = pls[key][row]
            with np.errstate(divide='ignore', invalid='ignore'):
                lr = np.log10(np.where(P > 0, P, 1e-30) / E)      # P_k <= 0 (backscatter): ratio below the scale
            lr = np.where(np.isfinite(lr), lr, np.nan)            # drawn at the -1 end, outlined by the P_k = 0 contour
            neg = negligible(P, E)
            nrm = Normalize(-1, 1)
            m = ax.pcolormesh(XP, ZG, np.where(neg, np.nan, lr), cmap='bwr', norm=nrm, shading='auto',
                              rasterized=True)
            if neg.any():
                ax.pcolormesh(XP, ZG, np.where(neg, 1.0, np.nan), cmap=ListedColormap(['#D9D9D9']),
                              shading='auto', rasterized=True)
                print(f'A panel {chr(97 + k)}: grey fraction {neg.mean():.3f}')
            if np.nanmin(P) < 0:
                ax.contour(XP, ZG, np.ma.masked_where(neg, P), levels=[0.0], colors='k', linewidths=0.6,
                           linestyles='--')
            mark_hole(ax, row)
            ax.text(0.995, 0.93, f'({chr(97 + k)}) {name}, {rl}', transform=ax.transAxes, ha='right', va='top',
                    bbox=dict(fc='w', ec='none', alpha=0.75, pad=1.0))
            ax.set_ylabel(r'$z/D$')
            ax.set_yticks([0, 1, 2])
            k += 1
    axs[-1].set_xlabel(r'$x/D$')
    fig.subplots_adjust(left=0.08, right=0.87, top=0.99, bottom=0.09, hspace=0.12)
    cax = fig.add_axes([0.89, 0.09, 0.018, 0.90])
    cb = fig.colorbar(m, cax=cax)
    cb.solids.set_alpha(1.0)
    cb.set_label(r'$\log_{10}(P_k/\varepsilon)$')
    return fig


def figure_B(plt, pls):
    from matplotlib.colors import LogNorm
    fig, axs = plt.subplots(6, 2, figsize=(7.2, 7.6), sharex=True)
    # layout: rows = plane x rung, columns would not fit 27 D; kept as 12 stacked full-width panels instead
    plt.close(fig)
    fig, axs = plt.subplots(12, 1, figsize=(7.2, 8.6), sharex=True)
    lo, hi = 1e-5, 1e-2
    k = 0
    for r, (row, rl) in enumerate(ROWS):
        for key, name in [('rans', 'RANS 50M'), ('ref', '193M 3DNS')]:
            for q, qn in enumerate((r'$P_k$', r'$\varepsilon$')):
                ax = axs[k]
                A = pls[key][row][q]
                m = ax.pcolormesh(XP, ZG, np.where(A > 0, A, lo), cmap='inferno', norm=LogNorm(lo, hi),
                                  shading='auto', rasterized=True)     # values <= 0 drawn at the scale floor
                if q == 0 and np.nanmin(A) < 0:
                    ax.contour(XP, ZG, A, levels=[0.0], colors=OKABE[2], linewidths=0.6, linestyles='--')
                mark_hole(ax, row)
                ax.text(0.995, 0.90, f'({chr(97 + k)}) {qn}, {name}, {rl}', transform=ax.transAxes, ha='right',
                        va='top', color='w', fontsize=plt.rcParams['font.size'] * 0.85)
                ax.set_yticks([0, 1])
                ax.set_ylabel(r'$z/D$')
                k += 1
    axs[-1].set_xlabel(r'$x/D$')
    fig.subplots_adjust(left=0.08, right=0.87, top=0.99, bottom=0.07, hspace=0.12)
    cax = fig.add_axes([0.89, 0.07, 0.018, 0.92])
    cb = fig.colorbar(m, cax=cax)
    cb.set_label(r'$(P_k,\ \varepsilon)\,D/(\rho_\infty U_\infty^3)$')
    return fig


def main():
    var = sys.argv[1] if len(sys.argv) > 1 else 'AB'
    if var == 'plan':
        run_ref(plan_only=True)
        return
    if var == 'rans':
        run_rans()
        return
    if var == 'ref':
        run_ref()
        return
    if var == 'check':
        check()
        return
    ref, rs = run_ref(), run_rans()
    pls = {'ref': planes(ref), 'rans': planes(rs)}
    print(f'threshold: |P| and eps both below {THRESH:g} rho U^3/D -> flat grey in variant A')
    for key in pls:
        for row, _ in ROWS:
            print(key, row, {k: round(v, 3) for k, v in medians(pls[key][row]).items()})
    plt = setup()
    if 'A' in var:
        fig = figure_A(plt, pls)
        fig.savefig(os.path.join(OUTDIR, 'fig_budget_maps.png'), dpi=300)
        print('saved fig_budget_maps.png')
    if 'B' in var:
        fig = figure_B(plt, pls)
        fig.savefig(os.path.join(OUTDIR, 'fig_budget_maps_B.png'), dpi=300)
        print('saved fig_budget_maps_B.png')


if __name__ == '__main__':
    main()
