#!/usr/bin/env python3
"""Time-batch delete-one jackknife of the mode C regional split, and dissipation-length check (2026-11-26).

Input: eps_points_c193b_w1_b010.npz (fc session point file; 193M batch b010 only). Nothing else is read.
Part 1: ln(mu_R,RANS/mu_R,193) and the 4-player exact Shapley split (rho, energy = k^1/2, L_al = mu_eff/(rho k^1/2), strain weights)
        per region, full record, then delete-one over the 10 flow-through batches (RANS fixed, 193M from jk_ fields):
        SE = sqrt((n-1)/n sum (theta_b - mean)^2). Mask: mode C (wall nodes excluded, z/D >= 0.3 on plate regions, far film x/D <= 20).
        The 5-player split (L_eff, alignment) needs per-batch L and cos, which the file does not hold: not computed.
Part 2: volume-weighted median / IQR of L_eps/D, L_al/D, L_eff/D (193M) and L_R/D (RANS), C = L_al/L_eps, fraction with eps <= 0.
Writes figures/code/data/meffR_se_time.txt.
usage: python3 -B estimator_timejk.py
"""
import itertools
import math

import numpy as np

F = '/mnt/Elements/DNS/analysis/asme2027_w1/eps_points_c193b_w1_b010.npz'
OUT = '/home/sk999/Documents/ASME_2027/figures/code/data/meffR_se_time.txt'
D = 0.00775
d = np.load(F)
names = [str(s) for s in d['region_names']]
reg = d['region']
xyz = d['xyz']
V = d['w'].astype(np.float64) / 8.0
g = lambda n: d[n].astype(np.float64)  # noqa: E731
rD, kD, LD, cD, SD, mD = g('rho193'), g('k193'), g('L193'), g('cos193'), g('S193'), g('mueff193')
rR, kR, mutR, SR = g('rhoR'), g('kR'), g('mutR'), g('SR')
eps = g('eps193_b010') / rD
jk = {n: d[n] for n in ('jk_mueff193', 'jk_S2_193', 'jk_k193', 'jk_rho193')}
REG = [('jet_0-3', 'near_jet_0_3'), ('film_3-10', 'wake_3_10'), ('film_10-25', 'wake_10_25')]
PLATE = {'near_jet_0_3', 'wake_3_10', 'wake_10_25'}

fin = np.all([np.isfinite(x) for x in (rD, kD, LD, cD, SD, mD, rR, kR, mutR, SR)], axis=0)
valid = fin & (kD > 0) & (kR > 0)          # wall nodes (k_193 = 0 or k_R <= 0) excluded
for nm in PLATE:
    valid &= ~((reg == names.index(nm)) & (xyz[:, 2] < 0.3 * D))
valid &= ~((reg == names.index('wake_10_25')) & (xyz[:, 0] > 20 * D))
print('file note:', str(d['note']))
print('xyz range x/D', xyz[:, 0].min() / D, xyz[:, 0].max() / D, ' z/D', xyz[:, 2].min() / D, xyz[:, 2].max() / D)

SUBS = list(itertools.product([False, True], repeat=4))   # rho, energy, L_al, strain weight ; False = RANS state? see below
# state True = 193M value, False = RANS value, as in estimator.py (value at subset s: players in s take the 193M side)? estimator.py:
# A[q][int(s[i])]: index 0 = 193M, 1 = RANS. s all False -> 193M; all True -> RANS. total = fv[all True] - fv[all False] = ln RANS/193.


def sums_for(idx, A):
    out = {}
    v = V[idx]
    for s in SUBS:
        x = [A[i][int(s[i])] for i in range(4)]
        den = v * x[3] ** 2
        num = den * x[0] * x[1] * x[2]
        out[s] = (num.sum(), den.sum())
    return out


def terms(sums):
    fv = {s: math.log(n / dn) for s, (n, dn) in sums.items()}
    tot = fv[(True,) * 4] - fv[(False,) * 4]
    phi = np.zeros(4)
    od = list(itertools.permutations(range(4)))
    for o in od:
        st = [False] * 4
        prev = fv[tuple(st)]
        for i in o:
            st[i] = True
            cur = fv[tuple(st)]
            phi[i] += (cur - prev) / len(od)
            prev = cur
    return np.concatenate([[tot], phi])   # total, rho, energy, L_al, S-wt


with np.errstate(all='ignore'):
    LalR = mutR / (rR * np.sqrt(kR))
    LalD = mD / (rD * np.sqrt(kD))
print('\nPART 1. mode C reproduction and time-batch jackknife')
print('columns: total rho energy L_al S-wt  (ln units)')
res = {}
for lab, nm in REG:
    ri = names.index(nm)
    idx = np.where(valid & (reg == ri))[0]
    RA = [rR[idx], np.sqrt(kR[idx]), LalR[idx], SR[idx]]
    full = terms(sums_for(idx, [(rD[idx], RA[0]), (np.sqrt(kD[idx]), RA[1]), (LalD[idx], RA[2]), (SD[idx], RA[3])]))
    # check: L_al from mueff equals L_eff * cos
    chk = np.nanmax(np.abs(LalD[idx] - LD[idx] * cD[idx]) / np.abs(LalD[idx]).clip(1e-30))
    t = []
    nbad = 0
    for b in range(10):
        r_, k_, m_, s2 = (jk[n][b][idx].astype(np.float64) for n in ('jk_rho193', 'jk_k193', 'jk_mueff193', 'jk_S2_193'))
        ok = np.isfinite(r_) & np.isfinite(k_) & np.isfinite(m_) & np.isfinite(s2) & (k_ > 0) & (r_ > 0)
        nbad += int((~ok).sum())
        r_, k_, m_, s2 = (np.where(ok, x, 1.0) for x in (r_, k_, m_, s2))
        wgt = np.where(ok, 1.0, 0.0)
        # points failing in a batch: set S2 = 0 (zero weight)
        s2 = s2 * wgt
        Sb = np.sqrt(np.maximum(s2, 0))
        Lb = m_ / (r_ * np.sqrt(k_))
        t.append(terms(sums_for(idx, [(r_, RA[0]), (np.sqrt(k_), RA[1]), (Lb, RA[2]), (Sb, RA[3])])))
    t = np.array(t)
    se = np.sqrt(9 / 10 * ((t - t.mean(0)) ** 2).sum(0))
    res[lab] = (full, se, t)
    print(f'{lab:11s} n={len(idx)} vol={V[idx].sum():.3e} m3  max|L_al - L_eff cos|/L_al = {chk:.1e}; jk points failing = {nbad}')
    print('   full   ' + ' '.join(f'{x:+7.3f}' for x in full) + f'   sum-terms-total {full[1:].sum() - full[0]:+.1e}')
    print('   jk mean' + ' '.join(f'{x:+7.3f}' for x in t.mean(0)))
    print('   SE     ' + ' '.join(f'{x:7.3f}' for x in se))
    print('   range  ' + ' '.join(f'[{a:+.3f},{b:+.3f}]' for a, b in zip(t.min(0), t.max(0))))
    print('   |value|/SE ' + ' '.join(f'{abs(a) / b:7.1f}' for a, b in zip(full, se)))

with open(OUT, 'w', encoding='utf-8') as f:
    f.write('# kind: time-batch jackknife\n')
    f.write('# SE of the total ln(mu_R,RANS / mu_R,193), delete-one over 10 flow-through batches (RANS fixed); source estimator_timejk.py\n')
    f.write('region SE density_SE energy_SE L_al_SE Swt_SE\n')
    for lab, _ in REG:
        se = res[lab][1]
        f.write(f'{lab} {se[0]:.3f} {se[1]:.3f} {se[2]:.3f} {se[3]:.3f} {se[4]:.3f}\n')
print('wrote', OUT)

# ------------------------------------------------------------------ Part 2
print('\nPART 2. dissipation length (same masks; volume weights dV; eps <= 0 points excluded from the L_eps statistics)')


def wq(x, w, q):
    o = np.argsort(x)
    c = np.cumsum(w[o])
    c /= c[-1]
    return np.interp(q, c, x[o])


def stat(x, w):
    return tuple(wq(x, w, q) for q in (0.25, 0.5, 0.75))


for lab, nm in REG:
    ri = names.index(nm)
    m = valid & (reg == ri)
    vt = V[m].sum()
    me = m & (eps > 0)
    fneg = V[m & ~(eps > 0)].sum() / vt
    with np.errstate(all='ignore'):
        Leps = kD ** 1.5 / eps
        C = LalD / Leps
    w = V[me]
    print(f'\n{lab}: volume {vt:.3e} m3, fraction eps <= 0 (excluded) {fneg:.4f}')
    print(f'   {"quantity":22s} {"q25":>9s} {"median":>9s} {"q75":>9s}   (/D)')
    for nmq, x, mm in (('L_eps 193M', Leps, me), ('L_al 193M', LalD, me), ('L_eff 193M', LD, me), ('L_R RANS', LalR, me),
                       ('L_al 193M (all m)', LalD, m), ('L_eff 193M (all m)', LD, m), ('L_R RANS (all m)', LalR, m),
                       ('0.09 L_eps', 0.09 * Leps, me)):
        q = stat(x[mm] / D, V[mm])
        print(f'   {nmq:22s} {q[0]:9.4f} {q[1]:9.4f} {q[2]:9.4f}')
    q = stat(C[me], V[me])
    print(f'   C = L_al/L_eps pointwise: q25 {q[0]:+.4f} median {q[1]:+.4f} q75 {q[2]:+.4f}')
    # integral forms
    Cint = (V[me] * LalD[me]).sum() / (V[me] * Leps[me]).sum()
    print(f'   C from volume integrals sum V L_al / sum V L_eps = {Cint:.4f}; ratio of medians L_al/L_eps = '
          f'{stat(LalD[me] / D, V[me])[1] / stat(Leps[me] / D, V[me])[1]:.4f}')
    print(f'   RANS L_R median / (0.09 L_eps median) = {stat(LalR[me] / D, V[me])[1] / (0.09 * stat(Leps[me] / D, V[me])[1]):.3f}; '
          f'L_R / L_al median = {stat(LalR[me] / D, V[me])[1] / stat(LalD[me] / D, V[me])[1]:.3f}')
    # fraction of volume with 193M L_al < 0 (counter-gradient)
    print(f'   volume fraction L_al <= 0 (193M): {V[m & (LalD <= 0)].sum() / vt:.3f}')
