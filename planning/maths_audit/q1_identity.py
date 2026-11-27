#!/usr/bin/env python3
"""Q1: is mu_eff = rho k^1/2 L_eff cos(theta) exact as coded?

(a) Code-faithful synthetic test: the operations of analysis/meff_c50b_c100b/dns_fields.py lines 56-88 are copied
    onto random Favre moments (rho, <rho u>, <rho u_i u_j>) and a random velocity-gradient tensor, and
    -a:S*/(2|S*|^2) is compared with rho k^1/2 L cos. Also: a = rho u''u'' - 2/3 rho k I from the raw sums equals
    the Favre definition; a:S* sign for a Boussinesq stress; |b| realizability bound.
(b) Cached 193M b010 points (meff_planes/part_b010_*.npz, 5.6 M points): rho k^1/2 L cos against the independently
    coded meff = -aS/(2 S2) (meff_planes_common.py line 149).
usage: python3 -B q1_identity.py
"""
import glob
import numpy as np

rng = np.random.default_rng(1)
N = 200000
# ---- (a) synthetic, code path of dns_fields.py ----
rho = rng.uniform(0.5, 1.5, N)
Uf = rng.normal(0, 30, (N, 3))
# Favre fluctuation tensor: random SPD  <u''_i u''_j>
A = rng.normal(0, 3, (N, 3, 3))
Rf = np.einsum('nij,nkj->nik', A, A)                         # u''u'' (SPD)
# raw stat3 sums: <rho u_i u_j> = rho (U_i U_j + u''_i u''_j), <rho u_i> = rho U_i
rhou = rho[:, None] * Uf
PAIRS = {(0, 0): 'rhouu', (1, 1): 'rhovv', (2, 2): 'rhoww', (0, 1): 'rhouv', (0, 2): 'rhouw', (1, 2): 'rhovw'}
m = {c: rho * (Uf[:, i] * Uf[:, j] + Rf[:, i, j]) for (i, j), c in PAIRS.items()}
U = rhou / rho[:, None]
a = np.empty((N, 3, 3))
for (i, j), c in PAIRS.items():
    a[:, i, j] = a[:, j, i] = m[c] - rho * U[:, i] * U[:, j]
k = 0.5 * np.trace(a, axis1=-2, axis2=-1) / rho
for i in range(3):
    a[:, i, i] -= 2.0 / 3.0 * rho * k
a_ref = rho[:, None, None] * (Rf - (np.trace(Rf, axis1=1, axis2=2) / 3)[:, None, None] * np.eye(3))
print('Favre a from raw sums vs rho(u\'\'u\'\' - 2/3 k I): max abs rel diff %.2e' % (np.abs(a - a_ref).max() / np.abs(a_ref).max()))
G = rng.normal(0, 1e3, (N, 3, 3))                            # dU_i/dx_j, compressible (trace != 0)
Sm = 0.5 * (G + np.swapaxes(G, -1, -2))
tr = np.trace(Sm, axis1=-2, axis2=-1)
for i in range(3):
    Sm[:, i, i] -= tr / 3
S2 = np.einsum('nij,nij->n', Sm, Sm)
aS = np.einsum('nij,nij->n', a, Sm)
an = np.sqrt(np.einsum('nij,nij->n', a, a))
cos = -aS / (an * np.sqrt(S2))
L = np.sqrt(k) * (an / (2 * rho * k)) / np.sqrt(S2)
meff = -aS / (2 * S2)
prod = rho * np.sqrt(k) * L * cos
print('(a) random tensors: max rel residual |rho k^.5 L cos - (-a:S*/(2|S*|^2))| / |.| = %.2e' % (np.abs(prod - meff) / np.abs(meff)).max())
# Boussinesq stress: a = -2 mu S*  ->  cos = 1, L = mu/(rho k^1/2), meff = mu
mu = rng.uniform(1e-3, 1e-1, N)
kk = rng.uniform(1, 50, N)
ab = -2 * mu[:, None, None] * Sm
anb = np.sqrt(np.einsum('nij,nij->n', ab, ab))
cb = -np.einsum('nij,nij->n', ab, Sm) / (anb * np.sqrt(S2))
Lb = np.sqrt(kk) * (anb / (2 * rho * kk)) / np.sqrt(S2)
print('    Boussinesq: cos-1 max %.1e ; L k^.5 rho/mu - 1 max %.1e' % (np.abs(cb - 1).max(), np.abs(Lb * np.sqrt(kk) * rho / mu - 1).max()))
# |a|^2 = |a_par|^2 + |a_perp|^2
ap = a + 2 * meff[:, None, None] * Sm
print('    a_perp:S* max |.|/(|a||S|) = %.1e ; |a|^2-|a_par|^2-|a_perp|^2 rel %.1e' % (
    (np.abs(np.einsum('nij,nij->n', ap, Sm)) / (an * np.sqrt(S2))).max(),
    np.abs(an**2 - (an * cos)**2 - np.einsum('nij,nij->n', ap, ap)).max() / (an**2).max()))
# ---- (b) cached 193M points ----
P = '/mnt/Elements/DNS/analysis/asme2027_figs/meff_planes/'
rels, bmax, nneg, ntot = [], 0.0, 0, 0
for f in sorted(glob.glob(P + 'part_b010_*.npz')):
    d = np.load(f)
    for s in range(4):
        r, kq, Lq, c, mq, Sq = [d[f'{q}{s}'].astype(np.float64) for q in ('rho', 'k', 'L', 'cos', 'meff', 'S')]
        ok = np.isfinite(mq) & (kq > 0) & (Sq > 0) & (mq != 0)
        if not ok.any():
            continue
        rels.append(np.abs(r * np.sqrt(kq) * Lq * c - mq)[ok] / np.abs(mq[ok]))
        bmax = max(bmax, float((Lq * Sq / np.sqrt(kq))[ok].max()))      # |b| = L |S*| / k^1/2
        nneg += int((c[ok] <= 0).sum()); ntot += int(ok.sum())
a_ = np.concatenate(rels)
print('(b) cached 193M b010 points: n = %d, max rel residual %.2e, median %.2e, 99.99 pct %.2e (float32 storage)' % (
    len(a_), a_.max(), np.median(a_), np.percentile(a_, 99.99)))
print('    max |b| = L|S*|/k^.5 = %.4f (realizability bound sqrt(2/3) = %.4f); cos<=0 fraction of these points %.3f' % (
    bmax, np.sqrt(2 / 3), nneg / ntot))
