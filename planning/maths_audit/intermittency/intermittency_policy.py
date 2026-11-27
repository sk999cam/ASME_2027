#!/usr/bin/env python3
"""Threshold policy numbers for the high-pass spanwise-stress intermittency detector (see intermittency.md).

Reads the caches of figures/code/_gamma_extract_v.py (z010 plane; pooled histograms of log10 d; fractions of time above 91
log10-threshold levels for T = 1, 2, 4 D/U; the mid-pitch line) and prints: Otsu thresholds and bimodality for two pooling
domains (plane x/D 0..20 all y; mid-pitch line x/D 0..20), mid-pitch gamma at x/D -2, 5, 10, 15, 20 for the policies i, ii, iii,
the two required checks, and the sensitivity to T and to the threshold. Output also written to intermittency_policy.txt.
Policies: i = each run its own plane-pooled Otsu threshold; ii = the plane-pooled threshold of the 193M for both runs;
iii = the mid-pitch-line-pooled threshold of the 193M for both runs.
usage: python3 intermittency_policy.py
"""
import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter1d

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, '/home/sk999/Documents/ASME_2027/figures/code')
import _gamma_extract_v as GV  # noqa: E402

TAGS = (('c193b', '193M'), ('c20bw', '20M LES'))
XQ = (-2, 5, 10, 15, 20)
out = []


def P(s=''):
    out.append(s)
    print(s, flush=True)


def moments(h, edges):
    c = 0.5 * (edges[1:] + edges[:-1])
    p = h / h.sum()
    m = (p * c).sum()
    s2 = (p * (c - m) ** 2).sum()
    sk = (p * (c - m) ** 3).sum() / s2 ** 1.5
    ku = (p * (c - m) ** 4).sum() / s2 ** 2 - 3.0
    n = h.sum()
    return (sk ** 2 + 1) / (ku + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3))), sk, ku


def valley(h, edges, rebin=5, sigma=2, sep=0.5):
    """valley / lower peak of the log10 d histogram (rebinned to 0.05 decade, Gaussian-smoothed over `sigma` bins); the two
    highest local maxima at least `sep` decades apart. nan if the smoothed histogram has one maximum only."""
    n = len(h) // rebin * rebin
    hr = h[:n].reshape(-1, rebin).sum(1).astype(float)
    c = edges[:n + 1:rebin][:-1] + 0.5 * rebin * (edges[1] - edges[0])
    hs = gaussian_filter1d(hr, sigma)
    pk = sorted([i for i in range(1, len(hs) - 1) if hs[i] >= hs[i - 1] and hs[i] > hs[i + 1]], key=lambda i: -hs[i])
    if not pk:
        return np.nan, None
    oth = [i for i in pk[1:] if abs(c[i] - c[pk[0]]) > sep]
    if not oth:
        return np.nan, (round(float(c[pk[0]]), 2),)
    a, b = sorted((pk[0], oth[0]))
    return hs[a:b + 1].min() / min(hs[pk[0]], hs[oth[0]]), (round(float(c[a]), 2), round(float(c[b]), 2))


def line_gamma(g):
    """mid-pitch line y/D = +-3 (j = 300, j = 0 identical): g (NXR, NY) -> (NXR,)."""
    return 0.5 * (g[:, 0] + g[:, -1])


def at(g1, xq):
    return [float(np.interp(q, GV.XG, g1)) for q in xq]


def fmt(v):
    return ' '.join(f'{a:.3f}' for a in v)


Z = {t: GV.load(t, 'z010') for t, _ in TAGS}
TH = GV.thresholds()
TQ = int(np.argmin(np.abs(Z['c193b']['TS'] - 2.0)))
K = {'i': 'own plane', 'ii': '193M plane th', 'iii': '193M line th', 'band': '193M band |y| 2-3 th',
     'band1.5': '193M band |y| 1.5-3 th', 'band2.5': '193M band |y| 2.5-3 th'}
POLS = ('i', 'ii', 'iii', 'band')

P('1. Otsu on log10 d, T = 2 D/U, two pooling domains (x/D 0..20)')
P('domain, run | log10 th | th | eta (between-class variance / total variance) | low-class fraction | valley/low-peak | '
  'histogram peaks (log10 d) | Sarle BC | skewness | excess kurtosis')
for dom in ('plane', 'line', 'band2', 'band1.5', 'band2.5'):
    for t, nm in TAGS:
        r = TH['2'][f'{dom}_{t}']
        h = np.asarray(r['hist']) if dom != 'plane' else Z[t]['hist'][TQ]
        bc, sk, ku = moments(h, GV.EDGES)
        vl, pk = valley(h, GV.EDGES)
        P(f'{dom}, {nm} | {r["log10"]:.3f} | {10 ** r["log10"]:.3e} | {r["eta"]:.3f} | {r["w0"]:.3f} | {vl:.3f} | {pk} | '
          f'{bc:.3f} | {sk:.3f} | {ku:.3f}')
P(f'threshold ratio LES/193M: plane {10 ** (TH["2"]["plane_c20bw"]["log10"] - TH["2"]["plane_c193b"]["log10"]):.3f}, '
  f'line {10 ** (TH["2"]["line_c20bw"]["log10"] - TH["2"]["line_c193b"]["log10"]):.3f}')
P('Sarle BC > 0.555 indicates bimodality (uniform distribution 0.555).')

G = {}
P('\n2. mid-pitch gamma(x), line y/D = +-3, T = 2 D/U; x/D ' + ' '.join(str(q) for q in XQ))
for pol in POLS:
    for t, nm in TAGS:
        lt = GV.policy_logth(t, 2.0, pol)
        g = GV.gamma_at(Z[t]['frac'][TQ], lt)
        G[(pol, t)] = g
        P(f'policy {pol} ({K[pol]}), {nm} (log10 th {lt:.3f}): ' + fmt(at(line_gamma(g), XQ)))
jb = np.where(np.abs(GV.YG) >= 2.81 - 1e-9)[0]
P('band |y/D| >= 2.81 mean (as fig_eta_x):')
for pol in POLS:
    for t, nm in TAGS:
        P(f'policy {pol}, {nm}: ' + fmt(at(G[(pol, t)][:, jb].mean(1), XQ)))

P('\n3. required checks (x/D -3..-1.5, all y: gamma < 0.02; |y/D| < 1, x/D 5..20: gamma > 0.95)')
xa = (GV.XG >= -3.0 - 1e-9) & (GV.XG <= -1.5)
xc = (GV.XG >= 5.0 - 1e-9) & (GV.XG <= 20.0 + 1e-9)
yc = np.abs(GV.YG) < 1.0
for pol in POLS:
    for t, nm in TAGS:
        g = G[(pol, t)]
        a = g[xa]
        c = g[np.ix_(xc, yc)]
        P(f'policy {pol}, {nm}: approach mean {a.mean():.4f}, p99 {np.percentile(a, 99):.4f}, max {a.max():.4f}; '
          f'core mean {c.mean():.3f}, min {c.min():.3f}, p1 {np.percentile(c, 1):.3f}, fraction of points > 0.95 '
          f'{(c > 0.95).mean():.3f}')
P('core gamma(x) along y/D = 0 and the |y/D| < 1 mean, policy ii vs iii, x/D 5 10 15 20:')
for pol in ('iii', 'band'):
    for t, nm in TAGS:
        g = G[(pol, t)]
        P(f'policy {pol}, {nm}: y=0 ' + fmt(at(g[:, np.argmin(np.abs(GV.YG))], (5, 10, 15, 20))) + ' | |y|<1 mean '
          + fmt(at(g[:, yc].mean(1), (5, 10, 15, 20))))

P('\n4. sensitivity to T (running-mean length); threshold of the same T from the policy named')
for pol in ('iii', 'band'):
    for q, T in enumerate(Z['c193b']['TS']):
        for t, nm in TAGS:
            lt = GV.policy_logth(t, float(T), pol)
            g = line_gamma(GV.gamma_at(Z[t]['frac'][q], lt))
            c = GV.gamma_at(Z[t]['frac'][q], lt)[np.ix_(xc, yc)]
            P(f'policy {pol}, T = {T:g} D/U, {nm}: log10 th {lt:.3f}; ' + fmt(at(g, XQ)) + f' | core mean {c.mean():.3f} min {c.min():.3f}')

P('\n5. sensitivity to the threshold (T = 2 D/U, factor on the policy-iii and policy-ii threshold, applied to both runs)')
for pol in ('band', 'iii'):
    for fac in (0.5, 1.0, 2.0):
        for t, nm in TAGS:
            lt = GV.policy_logth(t, 2.0, pol) + np.log10(fac)
            g = GV.gamma_at(Z[t]['frac'][TQ], lt)
            c = g[np.ix_(xc, yc)]
            P(f'policy {pol} x{fac:g} (log10 th {lt:.3f}), {nm}: ' + fmt(at(line_gamma(g), XQ))
              + f' | approach max {g[xa].max():.4f} | core mean {c.mean():.3f} min {c.min():.3f}')

P('\n6. first x/D (0..20) at which the mid-pitch gamma exceeds 0.5, T = 2')
for pol in POLS:
    for t, nm in TAGS:
        gl = line_gamma(G[(pol, t)])
        m = (GV.XG >= 0) & (GV.XG <= 20.0 + 1e-9) & (gl > 0.5)
        P(f'policy {pol}, {nm}: ' + (f'{GV.XG[m][0]:.2f}' if m.any() else 'not reached'))
P('\n7. sensitivity to the band limits (T = 2 D/U; 193M band threshold applied to both runs): x/D 5 10 15 20')
for pol in ('band1.5', 'band', 'band2.5'):
    for t, nm in TAGS:
        lt = GV.policy_logth(t, 2.0, pol)
        g = GV.gamma_at(Z[t]['frac'][TQ], lt)
        c = g[np.ix_(xc, yc)]
        P(f'{pol} (log10 th {lt:.3f}), {nm}: ' + fmt(at(line_gamma(g), (5, 10, 15, 20)))
          + f' | approach mean {g[xa].mean():.4f} max {g[xa].max():.4f} | core mean {c.mean():.3f} min {c.min():.3f}')
open(os.path.join(HERE, 'intermittency_policy.txt'), 'w').write('\n'.join(out) + '\n')
