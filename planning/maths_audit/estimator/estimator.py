#!/usr/bin/env python3
"""Regional gradient-weighted least-squares (LS) estimator of the effective eddy viscosity and its exact splits.

Input: pts.npz (scratchpad, written by ../audit_extract.py). Read-only on everything else. Output on stdout
(saved as estimator.log). No new 193M read.

Per point  mu_eff = -a:S*/(2|S*|^2) = rho k^1/2 L_eff cos ;  RANS: cos = 1, mu_t = rho k^1/2 L_R.
Regional   mu_R = sum V (-a:S*) / sum V 2|S*|^2 = sum V |S*|^2 mu_eff / sum V |S*|^2 .

Checks (section A), splits (section B: weighted geometric mean, Shapley), table with block jackknife (section C).
"""
import itertools
import math
import sys

import numpy as np

d = np.load(sys.argv[1])
names = [str(s) for s in d['regnames']]
LAB = {'hole_metering': 'hole, metering', 'hole_diffuser': 'hole, diffuser', 'exit': 'hole exit', 'near_jet_0_3': 'jet 0-3',
       'wake_3_10': 'film 3-10', 'wake_10_25': 'film 10-25', 'upstream_bl': 'approach BL'}
ORDER = ['near_jet_0_3', 'wake_3_10', 'wake_10_25', 'exit', 'hole_diffuser', 'hole_metering', 'upstream_bl']
SMIN = 0.05 * 86.2812 / 0.00775
reg, blk = d['reg'], d['blk']
V = d['w'].astype(np.float64) / 8.0
g = lambda n: d[n].astype(np.float64)  # noqa: E731
rD, kD, LD, cD, SD = g('c193_rho'), g('c193_k'), g('c193_L'), g('c193_cos'), g('c193_S')
rR, kR, LR, mut, SR = g('rans_rho'), g('rans_k'), g('rans_L'), g('rans_mut'), g('S_r')
fin = np.all([np.isfinite(x) for x in (rD, kD, LD, cD, SD, rR, kR, LR, mut, SR)], axis=0)
# modes: A = finite, k_193 >= 0, k_R > 0 (previous run); B = A without wall nodes (k_193 = 0 or k_R <= 0; cos set to 0 there);
#        C = B and, on plate regions (jet, film 3-10, film 10-25, approach layer), z/D >= 0.3 (E10 cut). Hole regions: no wall-distance cut yet.
MODE = sys.argv[2] if len(sys.argv) > 2 else 'A'
XYZ = np.load(sys.argv[3])['xyz'] if len(sys.argv) > 3 else None
PLATE = [names.index(r) for r in ('near_jet_0_3', 'wake_3_10', 'wake_10_25', 'upstream_bl')]
base = fin & (kD >= 0) & (kR > 0)
wall = base & ~((kD > 0) & (kR > 0))
zcut = np.zeros(len(reg), bool)
if MODE == 'C':
    assert XYZ is not None
    zcut = np.isin(reg, PLATE) & (XYZ[:, 2] < 0.3)
valid = base.copy()
if MODE in ('B', 'C'):
    valid &= ~wall
if MODE == 'C':
    valid &= ~zcut
# optional 4th argument: x/D limit of the far-film region (owner rule 2026-11-26: nothing beyond x/D 20). Default none = unchanged.
XLIM = float(sys.argv[4]) if len(sys.argv) > 4 else None
if XLIM is not None:
    assert XYZ is not None
    valid &= ~((reg == names.index('wake_10_25')) & (XYZ[:, 0] > XLIM))
    print('far film (region wake_10_25) restricted to x/D <=', XLIM)
print('MODE', MODE, ' valid fraction of all points', float(valid.mean()))
with np.errstate(all='ignore'):
    muD = rD * np.sqrt(np.maximum(kD, 0)) * LD * cD
    muR = mut
print('RANS identity mu_t = rho k^1/2 L_R: max rel. err', float(np.nanmax(np.abs(rR * np.sqrt(kR) * LR / mut - 1)[valid])))


def hdr(t):
    print('\n' + '=' * 110 + '\n' + t + '\n' + '=' * 110)


def LS(mu, S, m):
    return float((V[m] * S[m] ** 2 * mu[m]).sum() / (V[m] * S[m] ** 2).sum())


# ------------------------------------------------------------------ A. maths checks
hdr('A. Checks of the LS properties on the 193M points (scalar form; tensor form on random tensors below)')
rng = np.random.default_rng(1)
# tensor test: minimiser of sum V |a + 2 mu S|^2 is mu_R ; production identity
n = 20000
A = rng.normal(size=(n, 3, 3)); A = 0.5 * (A + A.transpose(0, 2, 1)); A -= np.trace(A, axis1=1, axis2=2)[:, None, None] / 3 * np.eye(3)
G = rng.normal(size=(n, 3, 3)); S = 0.5 * (G + G.transpose(0, 2, 1)); S -= np.trace(S, axis1=1, axis2=2)[:, None, None] / 3 * np.eye(3)
Vv = rng.uniform(0.5, 2, n)
aS = np.einsum('nij,nij->n', A, S); S2 = np.einsum('nij,nij->n', S, S)
mu_ls = (Vv * -aS).sum() / (2 * (Vv * S2).sum())
obj = lambda mu: float((Vv * np.einsum('nij,nij->n', A + 2 * mu * S, A + 2 * mu * S)).sum())  # noqa: E731
print(f'1a random tensors: mu_R = {mu_ls:+.6e}; objective at mu_R*(1+-1e-3, 1e-2, 1e-1) minus minimum = ' +
      ', '.join(f'{obj(mu_ls * (1 + e)) - obj(mu_ls):.3e}' for e in (-0.1, -0.01, -0.001, 0.001, 0.01, 0.1)))
w_ = Vv * S2
mu_i = -aS / (2 * S2)
print(f'1b random tensors: sum V(-a:S) = {(Vv * -aS).sum():.9e}  sum V 2 mu_R |S|^2 = {(Vv * 2 * mu_ls * S2).sum():.9e}')
print(f'1c random tensors: sum w mu_i / sum w = {(w_ * mu_i).sum() / w_.sum():+.9e} against mu_R {mu_ls:+.9e}')
# point data, 193M: |a| = 2 rho k^1/2 L_eff |S|  ; p = -a:S = 2 mu_eff S^2 ; objective per point |a|^2 - 4 mu p + 4 mu^2 S^2
print('\n193M points: LS minimum, explained fraction, identity residuals (mask: finite, k >= 0, k_R > 0)')
print(f'{"region":15s} {"mu_R [Pa s]":>12s} {"obj(mu_R)":>11s} {"min of 9-pt scan":>16s} {"scan arg/mu_R":>13s} {"prod. resid":>11s} {"sum(w mu)/sum(w)":>16s} {"explained |a|^2":>15s}')
for rn in ORDER:
    ri = names.index(rn)
    m = valid & (reg == ri)
    a2 = (2 * rD * np.sqrt(np.maximum(kD, 0)) * LD * SD) ** 2
    p = 2 * muD * SD ** 2
    mu = LS(muD, SD, m)
    f = lambda x: float((V[m] * (a2[m] - 4 * x * p[m] + 4 * x ** 2 * SD[m] ** 2)).sum())  # noqa: E731
    scan = [mu * (1 + e) for e in np.linspace(-0.2, 0.2, 9)]
    vals = [f(x) for x in scan]
    wsum = (V[m] * SD[m] ** 2)
    prod_res = abs((V[m] * p[m]).sum() - (V[m] * 2 * mu * SD[m] ** 2).sum()) / abs((V[m] * p[m]).sum())
    expl = ((V[m] * p[m]).sum()) ** 2 / ((V[m] * SD[m] ** 2).sum() * (V[m] * a2[m]).sum())
    print(f'{LAB[rn]:15s} {mu:12.4e} {f(mu):11.4e} {min(vals):16.4e} {scan[int(np.argmin(vals))] / mu:13.3f} {prod_res:11.1e} '
          f'{(wsum * muD[m]).sum() / wsum.sum():16.4e} {expl:15.3f}')

hdr('A2. LS against thresholded mean of ln (current): which points carry the weight')
print(f'{"region":15s} {"LS D":>10s} {"<ln mu> D":>10s} {"LS R":>10s} {"<ln mu> R":>10s} {"ln LS ratio":>11s} {"<ln> ratio":>10s}  '
      f'{"w share cos<=0":>14s} {"w share S<Smin":>14s} {"w share cos<.3":>14s} {"w share k=0":>11s}')
for rn in ORDER:
    ri = names.index(rn)
    m = valid & (reg == ri)
    mA = m & (muD > 0) & (cD > 0) & (SD > SMIN) & (muR > 0)
    wD = V * SD ** 2
    lnm = lambda mu, mm: float((V[mm] * np.log(mu[mm])).sum() / V[mm].sum())  # noqa: E731
    sh = lambda c: float(wD[m & c].sum() / wD[m].sum())  # noqa: E731
    print(f'{LAB[rn]:15s} {LS(muD, SD, m):10.3e} {np.exp(lnm(muD, mA)):10.3e} {LS(muR, SR, m):10.3e} {np.exp(lnm(muR, mA)):10.3e} '
          f'{np.log(LS(muR, SR, m) / LS(muD, SD, m)):+11.3f} {lnm(muR, mA) - lnm(muD, mA):+10.3f}  '
          f'{sh(cD <= 0):14.3f} {sh(SD <= SMIN):14.3f} {sh((cD < 0.3)):14.3f} {sh(kD == 0):11.3f}')

# ------------------------------------------------------------------ B. splits
FACT = ['rho', 'k', 'Leff', 'cos', 'S']   # players


SUBS = list(itertools.product([False, True], repeat=5))
blocks = np.unique(blk)


def region_sums(ri, extra=None, Dd=None):
    """per block and per state: numerator, denominator over valid points of region ri."""
    m = valid & (reg == ri)
    if extra is not None:
        m = m & extra
    idx = np.where(m)[0]
    v, bb = V[idx], blk[idx]
    Dd = Dd or dict(rho=rD, k=kD, L=LD, c=cD, S=SD)      # 193M-side arrays; per-batch (delete-one) fields can be passed here
    A = dict(rho=(Dd['rho'][idx], rR[idx]), k=(np.sqrt(np.maximum(Dd['k'][idx], 0)), np.sqrt(np.maximum(kR[idx], 0))), L=(Dd['L'][idx], LR[idx]),
             c=(Dd['c'][idx], np.ones(len(idx))), S=(Dd['S'][idx], SR[idx]))
    keys = ['rho', 'k', 'L', 'c', 'S']
    out = {b: {} for b in np.unique(bb)}
    for s in SUBS:
        x = {q: A[q][int(s[i])] for i, q in enumerate(keys)}
        den = v * x['S'] ** 2
        num = den * x['rho'] * x['k'] * x['L'] * x['c']
        for b in out:
            sel = bb == b
            out[b][s] = (num[sel].sum(), den[sel].sum())
    return out


def shapley(fv, players):
    """fv: dict subset(tuple of 5 bools)->ln mu_R. players: list of lists of indices (groups move together)."""
    npl = len(players)
    phi = np.zeros(npl)
    rng_ = {i: [np.inf, -np.inf] for i in range(npl)}
    # exact Shapley by enumerating orders
    orders = list(itertools.permutations(range(npl)))
    for o in orders:
        st = [False] * 5
        prev = fv[tuple(st)]
        for i in o:
            for j in players[i]:
                st[j] = True
            cur = fv[tuple(st)]
            dlt = cur - prev
            phi[i] += dlt / len(orders)
            rng_[i][0] = min(rng_[i][0], dlt)
            rng_[i][1] = max(rng_[i][1], dlt)
            prev = cur
    return phi, rng_


def lnmu(sums, blocks_use, s):
    N = sum(sums[b][s][0] for b in blocks_use)
    Dn = sum(sums[b][s][1] for b in blocks_use)
    return math.log(N / Dn) if N / Dn > 0 else np.nan


P5 = [[0], [1], [2], [3], [4]]
P4 = [[0], [1], [2, 3], [4]]     # L_al = L_eff cos as one player
hdr('B. Exact additive split of ln(mu_R,RANS / mu_R,193), Shapley over {rho, k^1/2, L_eff, cos, strain weights}; value = ln of regional LS viscosity')
RES = {}
for rn in ORDER:
    ri = names.index(rn)
    sums = region_sums(ri)
    bl = list(sums)

    def terms(bu):
        fv = {s: lnmu(sums, bu, s) for s in SUBS}
        tot = fv[(True,) * 5] - fv[(False,) * 5]
        p5, r5 = shapley(fv, P5)
        p4, r4 = shapley(fv, P4)
        # sequential order rho -> k -> L_eff -> cos -> S
        seq = []
        st = [False] * 5
        prev = fv[tuple(st)]
        for i in range(5):
            st[i] = True
            seq.append(fv[tuple(st)] - prev)
            prev = fv[tuple(st)]
        return dict(tot=tot, p5=p5, p4=p4, r5=r5, r4=r4, seq=np.array(seq))
    full = terms(bl)
    # delete-one jackknife over spatial blocks
    jk = [terms([b for b in bl if b != b0]) for b0 in bl] if len(bl) > 1 else []
    RES[rn] = (full, jk, bl)

print('All terms are ln units. energy = 1/2 ln k part. S-weights term = effect of weighting the RANS by its own |S*|^2 instead of the 193M one.')
print(f'{"region":15s} {"nblk":>4s} {"total":>7s} | 5-player Shapley: {"rho":>6s} {"energy":>6s} {"L_eff":>6s} {"align":>6s} {"S-wt":>6s} {"L_al=L+al":>9s} {"sum-tot":>8s} | 4-player: {"rho":>6s} {"energy":>6s} {"L_al":>6s} {"S-wt":>6s}')
for rn in ORDER:
    f, jk, bl = RES[rn]
    p5, p4 = f['p5'], f['p4']
    print(f'{LAB[rn]:15s} {len(bl):4d} {f["tot"]:+7.3f} |{"":18s}{p5[0]:+6.3f} {p5[1]:+6.3f} {p5[2]:+6.3f} {p5[3]:+6.3f} {p5[4]:+6.3f} {p5[2] + p5[3]:+9.3f} {p5.sum() - f["tot"]:+8.0e} |{"":10s} '
          f'{p4[0]:+6.3f} {p4[1]:+6.3f} {p4[2]:+6.3f} {p4[3]:+6.3f}')
print('\n5 players: min/max over the 120 orders of the marginal contribution of each player (order dependence), and the sequential order rho,k,L_eff,cos,S')
print(f'{"region":15s} ' + ' '.join(f'{x:>22s}' for x in ('rho', 'energy', 'L_eff', 'align', 'S-wt')) + '   sequential (rho,k,L,cos,S)')
for rn in ORDER:
    f, jk, bl = RES[rn]
    print(f'{LAB[rn]:15s} ' + ' '.join(f'[{f["r5"][i][0]:+6.3f},{f["r5"][i][1]:+6.3f}]'.rjust(22) for i in range(5)) + '   ' +
          ' '.join(f'{x:+.3f}' for x in f['seq']))
print('\n4 players (L_al as one player, 24 orders): min/max')
for rn in ORDER:
    f, jk, bl = RES[rn]
    print(f'{LAB[rn]:15s} ' + ' '.join(f'[{f["r4"][i][0]:+6.3f},{f["r4"][i][1]:+6.3f}]'.rjust(22) for i in range(4)) + '   (rho, energy, L_al, S-wt)')

hdr('B(i). Weighted geometric means, weights w = V |S*|^2 (193M strain, common to both sides; per-point exactly additive). Mask mu_eff > 0 (cos > 0) needed.')
print(f'{"region":15s} {"ln LS ratio":>11s} {"sum terms":>9s} {"rho":>6s} {"energy":>6s} {"L_eff":>6s} {"align":>6s} {"L_al":>6s} {"gap":>6s} | own-weights: {"ln GM ratio":>11s} {"gap":>6s} | volume of valid kept (mask) {"":s}')
GM = {}
for rn in ORDER:
    ri = names.index(rn)
    m = valid & (reg == ri)
    mk = m & (muD > 0) & (cD > 0) & (kD > 0) & (LD > 0)
    w = V * SD ** 2
    wm = lambda x: float((w[mk] * x[mk]).sum() / w[mk].sum())  # noqa: E731
    with np.errstate(all='ignore'):
        trho, ten, tL, tal = wm(np.log(rR / rD)), wm(0.5 * np.log(kR / kD)), wm(np.log(LR / LD)), wm(-np.log(cD))
    tot = trho + ten + tL + tal
    lsr = np.log(LS(muR, SR, m) / LS(muD, SD, m))
    # own weights
    wR, wD_ = V * SR ** 2, V * SD ** 2
    gmR = float((wR[mk] * np.log(muR[mk])).sum() / wR[mk].sum())
    gmD = float((wD_[mk] * np.log(muD[mk])).sum() / wD_[mk].sum())
    kept = float(w[mk].sum() / w[m].sum())
    GM[rn] = dict(lsr=lsr, tot=tot, rho=trho, en=ten, L=tL, al=tal, gap=lsr - tot, kept=kept, gm_own=gmR - gmD)
    print(f'{LAB[rn]:15s} {lsr:+11.3f} {tot:+9.3f} {trho:+6.3f} {ten:+6.3f} {tL:+6.3f} {tal:+6.3f} {tL + tal:+6.3f} {lsr - tot:+6.3f} | {"":12s}{gmR - gmD:+11.3f} {lsr - (gmR - gmD):+6.3f} | strain-weight share kept {kept:.3f}')

# ------------------------------------------------------------------ C. table
hdr('C. Recommended estimator per region: LS value ratio with exact Shapley split (5 players, then L_al); block delete-one jackknife SE')
print('Jackknife = delete-one over the 22 spatial blocks of the point set (blocks that hold points of the region), NOT over time batches. The delete-one-flow-through fields were not read.')
print('SE_jk = sqrt((n-1)/n sum (t_i - mean)^2). Blocks are unequal in size and neighbouring blocks are correlated; read it as a sensitivity to spatial sub-sampling.')
print('Excluded volume = share of region volume with a non-finite value on either side (50M buffer / outside 193M) or k_193 < 0; k_193 = 0 wall nodes stay in (mu = 0, weight |S*|^2 retained).')
print(f'{"region":15s} {"mu_R,193":>9s} {"mu_R,RANS":>9s} {"ln ratio":>8s} {"energy":>14s} {"rho":>14s} {"L_eff":>14s} {"align":>14s} {"L_al":>14s} {"S-wt":>14s}   {"nonfinite":>9s} {"wall":>9s} {"z cut":>9s} {"cos<=0 vol":>10s} {"k=0 vol":>8s}')
rows = {}
for rn in ORDER:
    ri = names.index(rn)
    f, jk, bl = RES[rn]
    def se(fn):
        v = np.array([fn(j) for j in jk]) if jk else np.array([np.nan])
        nn = len(v)
        return math.sqrt((nn - 1) / nn * ((v - v.mean()) ** 2).sum()) if nn > 1 else float('nan')
    m = valid & (reg == ri)
    allr = reg == ri
    ex = V[allr & ~fin | allr & fin & (kD < 0)].sum() / V[allr].sum()
    exw = V[allr & wall].sum() / V[allr].sum() if MODE in ('B', 'C') else 0.0
    exz = V[allr & zcut & ~wall].sum() / V[allr].sum()
    cneg = V[m & (cD <= 0)].sum() / V[m].sum()
    k0 = V[m & (kD == 0)].sum() / V[m].sum()
    cell = lambda val, fn: f'{val:+6.3f}+-{se(fn):.3f}'  # noqa: E731
    row = [LAB[rn], LS(muD, SD, m), LS(muR, SR, m), f['tot'], f['p5'][1], f['p5'][0], f['p5'][2], f['p5'][3], f['p5'][2] + f['p5'][3], f['p5'][4]]
    print(f'{LAB[rn]:15s} {LS(muD, SD, m):9.3e} {LS(muR, SR, m):9.3e} {cell(f["tot"], lambda j: j["tot"]):>14s} {cell(f["p5"][1], lambda j: j["p5"][1]):>14s} '
          f'{cell(f["p5"][0], lambda j: j["p5"][0]):>14s} {cell(f["p5"][2], lambda j: j["p5"][2]):>14s} {cell(f["p5"][3], lambda j: j["p5"][3]):>14s} '
          f'{cell(f["p5"][2] + f["p5"][3], lambda j: j["p5"][2] + j["p5"][3]):>14s} {cell(f["p5"][4], lambda j: j["p5"][4]):>14s}   {ex:9.3f} {exw:9.3f} {exz:9.3f} {cneg:10.3f} {k0:8.4f}')
print('\nRANS cos = 1 everywhere (checked: min = max = 1). The alignment term is then the full -ln of the 193M strain-weighted alignment; it is not -<ln cos> of the old mean-of-ln.')


# ------------------------------------------------------------------ D. sensitivity and hierarchical variant
hdr('D. Wall nodes (k_193 = 0, mu_eff = 0 by construction) and hierarchical split')
print('Strain-weight share (V |S*|^2_193 share of the region) held by k_193 = 0 nodes, and their volume share:')
for rn in ORDER:
    ri = names.index(rn)
    m = valid & (reg == ri)
    w = V * SD ** 2
    print(f'  {LAB[rn]:15s} weight share {w[m & (kD == 0)].sum() / w[m].sum():.3f}   volume share {V[m & (kD == 0)].sum() / V[m].sum():.4f}   '
          f'RANS weight share (V S_R^2) at the same nodes {(V * SR ** 2)[m & (kD == 0)].sum() / (V * SR ** 2)[m].sum():.3f}')
print('\nHierarchical split (6 numbers): ln ratio = [ln mu~_R,RANS - ln mu_R,193] (physical factors on the 193M strain weights; 4-player Shapley) + [ln mu_R,RANS - ln mu~_R,RANS] (weight redistribution).')
print('  mu~ = RANS regional value with the 193M weights V|S*|_193^2. Columns: all valid points | without k_193 = 0 nodes.')
print(f'{"region":15s} | {"total":>7s} {"energy":>7s} {"rho":>7s} {"L_eff":>7s} {"align":>7s} {"L_al":>7s} {"S-wt":>7s} | {"total":>7s} {"energy":>7s} {"rho":>7s} {"L_eff":>7s} {"align":>7s} {"L_al":>7s} {"S-wt":>7s}   Shapley5 S-wt (no k0)')
for rn in ORDER:
    ri = names.index(rn)
    row = f'{LAB[rn]:15s} |'
    extra5 = ''
    for extra in (None, kD > 0):
        sums = region_sums(ri, extra)
        bl = list(sums)
        fv = {s: lnmu(sums, bl, s) for s in SUBS}
        tot = fv[(True,) * 5] - fv[(False,) * 5]
        f4 = {s: fv[tuple(list(s[:4]) + [False])] for s in SUBS}
        pl = [[0], [1], [2], [3]]
        phi = np.zeros(4)
        ords = list(itertools.permutations(range(4)))
        for o in ords:
            st = [False] * 5
            prev = f4[tuple(st)]
            for i in o:
                st[i] = True
                cur = f4[tuple(st)]
                phi[i] += (cur - prev) / len(ords)
                prev = cur
        swt = tot - phi.sum()
        row += f' {tot:+7.3f} {phi[1]:+7.3f} {phi[0]:+7.3f} {phi[2]:+7.3f} {phi[3]:+7.3f} {phi[2] + phi[3]:+7.3f} {swt:+7.3f} |'
        p5, _ = shapley(fv, P5)
        extra5 = f'{p5[4]:+.3f}'
    print(row + f'   {extra5}')

hdr('E. Delete-one block range of the total ln(mu_R,RANS/mu_R,193), and the block holding the largest share of the 193M strain weight')
for rn in ORDER:
    ri = names.index(rn)
    f, jk, bl = RES[rn]
    t = np.array([j['tot'] for j in jk])
    m = valid & (reg == ri)
    wD_ = V * SD ** 2
    sh = {int(b): wD_[m & (blk == b)].sum() / wD_[m].sum() for b in bl}
    bmax = max(sh, key=sh.get)
    print(f'{LAB[rn]:15s} total {f["tot"]:+.3f}  delete-one range [{t.min():+.3f}, {t.max():+.3f}]  blocks {len(bl)}  largest weight share block {bmax}: {sh[bmax]:.2f}  '
          f'total without it {[j["tot"] for j, b in zip(jk, bl) if b == bmax][0]:+.3f}')

hdr('F. Headline 4-player Shapley (rho, energy, L_al = L_eff cos, strain weights), delete-one block jackknife SE')
print(f'{"region":15s} {"total":>14s} {"energy":>14s} {"rho":>14s} {"L_al":>14s} {"S-wt":>14s}  order range of energy / L_al (24 orders)')
for rn in ORDER:
    f, jk, bl = RES[rn]

    def se4(i):
        v = np.array([j['p4'][i] for j in jk])
        nn = len(v)
        return math.sqrt((nn - 1) / nn * ((v - v.mean()) ** 2).sum())
    st = np.array([j['tot'] for j in jk]); nn = len(st)
    set_ = math.sqrt((nn - 1) / nn * ((st - st.mean()) ** 2).sum())
    c = lambda val, s: f'{val:+6.3f}+-{s:.3f}'  # noqa: E731
    p = f['p4']
    print(f'{LAB[rn]:15s} {c(f["tot"], set_):>14s} {c(p[1], se4(1)):>14s} {c(p[0], se4(0)):>14s} {c(p[2], se4(2)):>14s} {c(p[3], se4(3)):>14s}  '
          f'[{f["r4"][1][0]:+.2f},{f["r4"][1][1]:+.2f}] / [{f["r4"][2][0]:+.2f},{f["r4"][2][1]:+.2f}]')


# ------------------------------------------------------------------ G. time-batch delete-one jackknife (optional)
# usage: DELETE_ONE_FILES="f1.npz f2.npz ..." python3 -B estimator.py pts.npz C xyz.npz
# Each file: delete-one (record minus batch i) 193M fields at the points of pts.npz, same order: c193_rho, c193_k, c193_L, c193_cos, c193_S, shape (N,).
# The point mask `valid` of the full record is kept fixed. SE_jk = sqrt((N-1)/N sum (t_i - mean)^2) over the N delete-one samples.
import os  # noqa: E402
if os.environ.get('DELETE_ONE_FILES'):
    hdr('G. Time-batch delete-one jackknife of the 4-player Shapley terms')
    files = os.environ['DELETE_ONE_FILES'].split()
    res = {rn: [] for rn in ORDER}
    for fn in files:
        z = np.load(fn)
        Dd = dict(rho=z['c193_rho'].astype(float), k=z['c193_k'].astype(float), L=z['c193_L'].astype(float), c=z['c193_cos'].astype(float), S=z['c193_S'].astype(float))
        for rn in ORDER:
            sums = region_sums(names.index(rn), None, Dd)
            bl = list(sums)
            fv = {s_: lnmu(sums, bl, s_) for s_ in SUBS}
            res[rn].append([fv[(True,) * 5] - fv[(False,) * 5]] + list(shapley(fv, P4)[0]))
    print(f'{"region":15s} {"total":>14s} {"rho":>14s} {"energy":>14s} {"L_al":>14s} {"S-wt":>14s}   (N = {len(files)} delete-one samples)')
    for rn in ORDER:
        a = np.array(res[rn]); nn = len(a)
        se = np.sqrt((nn - 1) / nn * ((a - a.mean(0)) ** 2).sum(0))
        full = RES[rn][0]
        vals = [full['tot']] + list(full['p4'])
        print(f'{LAB[rn]:15s} ' + ' '.join(f'{v:+6.3f}+-{e:.3f}'.rjust(14) for v, e in zip(vals, se)))
