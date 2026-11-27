#!/usr/bin/env python3
"""Compare the QCR form with STAR's modelled Reynolds stresses (cells_qcr.csv, Tank, read-only; byte-offset segments only).
usage: python3 -B validate_star_stress.py
For each segment: parse rows, locate them in the cached qcr gradient arrays (identity order of the export; checked on X,Y,Z,k),
a_STAR = rho (u'u'_ij - 2/3 k delta_ij) (columns 'Reynolds Stress' are J/kg = tau/rho, uu+vv+ww = 2k),
a_lin = -2 mu_t S*, comm = (O tau - tau O) with tau = 2 mu_t S*, O = 2W/|grad v|.
Fit a_STAR - a_lin = c comm (least squares over cells): c is the effective c_cr1 of the standard form; residual reported.
"""
import io
import numpy as np
import pandas as pd
import common as cm

F = '/mnt/Tank/STARCCM/film_cooling/rans_50M_Lm25b_BR15_FST_ZPG_p0dns_qcr/cells_qcr.csv'
cols = ['u', 'v', 'w', 'T', 'H', 'k', 'om', 'mut', 'g', 'uu', 'uv', 'uw', 'vv', 'vw', 'ww', 'rho', 'p', 'V', 'X', 'Y', 'Z']
import os
fs = os.path.getsize(F)
N = 53120567
Xc = np.load(f'{cm.C}/qcr/X.npy', mmap_mode='r')
tot = {'num': 0.0, 'den': 0.0, 'res': 0.0, 'lin': 0.0, 'n': 0}
rows_all = []
_i = np.arange(0, N, 20000)
_x, _z = np.asarray(Xc[_i]) / cm.D, np.asarray(np.load(f'{cm.C}/qcr/Z.npy', mmap_mode='r')[_i]) / cm.D
_film = _i[(_x >= 0) & (_x <= 20) & (_z >= 0.3) & (_z <= 2)]      # strided sample of the cache, locates film-dense row ranges
FRACS = (_film[np.linspace(0, len(_film) - 1, 12).astype(int)] + 3000) / N
for f in FRACS:
    off = int(f * fs)
    with open(F, 'rb') as fh:
        fh.seek(off)
        raw = fh.read(12_000_000)
    raw = raw[raw.index(b'\n') + 1: raw.rindex(b'\n')]
    df = pd.read_csv(io.BytesIO(raw), header=None, names=cols)
    x0 = df['X'].values[0]
    g = int(f * N)
    lo, hi = max(0, g - 4_000_000), min(N, g + 4_000_000)
    win = np.asarray(Xc[lo:hi], dtype=np.float64)
    cand = np.flatnonzero(np.abs(win - x0) < 1e-6 * max(abs(x0), 1e-3) + 1e-9) + lo
    start = None
    for c in cand:
        yy = np.asarray(np.load(f'{cm.C}/qcr/Y.npy', mmap_mode='r')[c:c + 3], dtype=np.float64)
        if np.allclose(yy, df['Y'].values[:3], rtol=1e-5, atol=1e-9):
            start = int(c)
            break
    if start is None:
        print('segment', f, 'row not located; skipped')
        continue
    n = len(df)
    idx = np.arange(start, start + n)
    d = cm.load('qcr', cm.NAMES, idx)
    kmatch = np.abs(d['k'] / df['k'].values - 1) < 1e-5            # row alignment per cell (k and |dX|/D)
    xmatch = np.abs(d['X'] - df['X'].values) / cm.D < 1e-5
    chk = 1.0 - float(np.mean(kmatch & xmatch))                       # fraction of misaligned rows
    G, S, Ss, W, gn = cm.tensors(d)
    z = d['Z'] / cm.D
    ok = kmatch & xmatch & (d['k'] > 0) & (z >= 0.3) & (z <= 2) & (d['X'] / cm.D >= 0) & (d['X'] / cm.D <= 20)
    if ok.sum() < 2000:
        print(f'segment f={f:.3f}: rows {start}..{start + n - 1}, only {ok.sum()} film cells; skipped', flush=True)
        continue
    aS = np.zeros((n, 3, 3))
    for (i, j), c in {(0, 0): 'uu', (0, 1): 'uv', (0, 2): 'uw', (1, 1): 'vv', (1, 2): 'vw', (2, 2): 'ww'}.items():
        aS[:, i, j] = aS[:, j, i] = df[c].values * df['rho'].values
    tr = np.trace(aS, axis1=1, axis2=2) / 3
    for i in range(3):
        aS[:, i, i] -= tr
    alin = -2 * d['mut'][:, None, None] * Ss
    comm = cm.qcr_part(d['mut'], Ss, W, gn, 1.0)        # c = 1 -> comm
    dlt = (aS - alin)[ok]
    cm_ = comm[ok]
    num = np.einsum('nij,nij->', dlt, cm_)
    den = np.einsum('nij,nij->', cm_, cm_)
    chat = num / den
    res = np.sqrt(np.einsum('nij,nij->', dlt - chat * cm_, dlt - chat * cm_) / np.einsum('nij,nij->', dlt, dlt))
    # per-cell c in cells with a sizeable QCR part
    pc = np.einsum('nij,nij->n', dlt, cm_) / np.maximum(np.einsum('nij,nij->n', cm_, cm_), 1e-300)
    cc = np.einsum('nij,nij->n', cm_, cm_)
    big = cc > np.quantile(cc, 0.5)
    print(f'segment f={f}: rows {start}..{start + n - 1}, misaligned row fraction (k, X) = {chk:.1e}, cells used {ok.sum()}, '
          f'c_cr1 fit = {chat:.4f}, rel. residual |dlt - c comm|/|dlt| = {res:.3f}, per-cell c median/p10/p90 = '
          f'{np.median(pc[big]):.4f}/{np.quantile(pc[big], .1):.4f}/{np.quantile(pc[big], .9):.4f}', flush=True)
    tot['num'] += num; tot['den'] += den; tot['n'] += int(ok.sum())
    # also: STAR a_qcr vs a_lin contraction; STAR stress trace check
    print('   (uu+vv+ww)/(2k) median', float(np.median((df['uu'] + df['vv'] + df['ww']).values[ok] / (2 * df['k'].values[ok]))),
          ' STAR a_perp:S* contraction / (|a_perp||S*|) max',
          float(np.max(np.abs(np.einsum('nij,nij->n', (aS - alin)[ok] - (np.einsum('nij,nij->n', aS[ok], Ss[ok]) / np.einsum('nij,nij->n', Ss[ok], Ss[ok]) + 2 * d['mut'][ok])[:, None, None] * Ss[ok], Ss[ok])))))
print('pooled c_cr1 fit', tot['num'] / tot['den'], 'cells', tot['n'])
