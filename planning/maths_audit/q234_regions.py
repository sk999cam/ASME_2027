#!/usr/bin/env python3
"""Q2, Q3, Q4 (and the region/point-set checks of Q6): regional terms of ln(mu_t/mu_eff^193) on the point set of
analysis/meff_rans50_c193b10 (every 2nd 50M point of the point-wise regions, 193M b010 interpolated there).

Input: pts.npz written by audit_extract.py (scratchpad). Output: text on stdout (saved as q234_regions.log).
No jackknife (the delete-one fields are 37 GB; not read). The mask is the one of errorbars.py lines 249-253 without
the delete-one condition, so f_V and values differ slightly from regions_se.txt.

Q2  length + alignment = <ln(L_RANS/L_al,193)>,  L_al = mu_eff/(rho k^1/2) = L_eff cos;  shared part <-ln cos_193>.
Q3  L_eps = C_mu k^3/2/eps (C_mu = 0.09), eps of the 50M 3DNS b006 budget (the 193M has no cached eps).
Q4  excluded volume, tail shares, mean of ln against ln of mean and the gradient-weighted LS regional mu_R.
"""
import sys
import numpy as np

CM = 0.09
D, U = 0.00775, 86.2812
SMIN = 0.05 * U / D
d = np.load(sys.argv[1])
names = [str(s) for s in d['regnames']]
LAB = {'hole_metering': 'hole, metering', 'hole_diffuser': 'hole, diffuser', 'exit': 'hole exit', 'near_jet_0_3': 'jet 0-3',
       'wake_3_10': 'film 3-10', 'wake_10_25': 'film 10-25', 'upstream_bl': 'approach BL'}
ORDER = ['near_jet_0_3', 'wake_3_10', 'wake_10_25', 'hole_metering', 'hole_diffuser', 'exit', 'upstream_bl']
reg, w = d['reg'], d['w'].astype(np.float64)
dV = w / 8.0
r1, k1, L1, c1, S1 = [d['c193_' + q].astype(np.float64) for q in ('rho', 'k', 'L', 'cos', 'S')]
rr, kr, Lr, mut = [d['rans_' + q].astype(np.float64) for q in ('rho', 'k', 'L', 'mut')]
Sr = d['S_r'].astype(np.float64)
mu1 = r1 * np.sqrt(np.maximum(k1, 0)) * L1 * c1
ok = (mut > 0) & (mu1 > 0) & (c1 > 0) & (S1 > SMIN) & (kr > 0)       # errorbars.py 249-250 (RANS cos = 1)
ok &= np.isfinite(mu1) & np.isfinite(mut) & np.isfinite(L1) & (k1 > 0)
with np.errstate(all='ignore'):
    ln = {'rho': np.log(rr / r1), 'energy': 0.5 * np.log(kr / k1), 'length': np.log(Lr / L1), 'align': -np.log(c1)}
    lnLal = np.log(Lr / (L1 * c1))                       # ln(L_RANS / L_al,193)
    mean = lambda x, m, ww=w: float((ww[m] * x[m]).sum() / ww[m].sum())  # noqa: E731


def hdr(t):
    print('\n' + '=' * 100 + '\n' + t + '\n' + '=' * 100)


hdr('Reproduction of regions_se.txt (no delete-one condition). Columns: total density energy length align f_V')
ref = {'film 3-10': (-0.501, 0.007, -0.351, -0.963, 0.806, 0.54), 'film 10-25': (-0.703, 0.004, -0.351, -0.928, 0.573, 0.48),
       'jet 0-3': (0.566, -0.002, 0.051, -0.126, 0.643, 0.53), 'hole exit': (-0.431, 0.028, -0.246, -0.641, 0.428, 0.98),
       'approach BL': (-0.499, 0.0, -0.967, -1.819, 2.287, 0.22), 'hole, diffuser': (-0.033, 0.004, -0.130, -0.680, 0.773, 0.87),
       'hole, metering': (0.704, -0.008, 0.533, -0.993, 1.173, 0.82)}
T = {}
for rn in ORDER:
    ri = names.index(rn)
    m = ok & (reg == ri)
    fv = w[m].sum() / w[reg == ri].sum()
    v = [mean(ln['rho'], m), mean(ln['energy'], m), mean(ln['length'], m), mean(ln['align'], m)]
    tot = sum(v)
    T[rn] = dict(rho=v[0], energy=v[1], length=v[2], align=v[3], total=tot, fV=fv, lal=mean(lnLal, m))
    print(f'{LAB[rn]:15s} mine: {tot:+.3f} {v[0]:+.3f} {v[1]:+.3f} {v[2]:+.3f} {v[3]:+.3f} {fv:.2f} | file: ' +
          ' '.join(f'{x:+.3f}' for x in ref[LAB[rn]][:5]) + f' {ref[LAB[rn]][5]:.2f}')

hdr('Q2: length and alignment are coupled. combined = length + align = <ln(L_RANS/L_al,193)>; shared = <ln 1/cos_193>')
print(f'{"region":15s} {"length":>8s} {"align":>8s} {"comb":>8s} {"<lnLal>":>8s} {"shared":>8s} {"energy":>8s} {"total":>8s}  '
      f'{"comb/total":>10s}  geo-mean cos193 = exp(-shared)')
for rn in ORDER:
    t = T[rn]
    comb = t['length'] + t['align']
    print(f'{LAB[rn]:15s} {t["length"]:+8.3f} {t["align"]:+8.3f} {comb:+8.3f} {t["lal"]:+8.3f} {t["align"]:+8.3f} {t["energy"]:+8.3f} '
          f'{t["total"]:+8.3f}  {comb / t["total"] if abs(t["total"]) > 0.05 else float("nan"):10.2f}  {np.exp(-t["align"]):.3f}')
print('length share of the shared ln(1/cos): |shared| / |length| = ' + ', '.join(
    f'{LAB[rn]} {abs(T[rn]["align"]) / abs(T[rn]["length"]):.2f}' for rn in ORDER))
# the 3-term split with L_al instead of L_eff: density, energy, L_al
print('\nTwo-factor alternative (no alignment term): total = density + energy + ln(L_RANS/L_al,193):')
for rn in ORDER:
    t = T[rn]
    print(f'  {LAB[rn]:15s} density {t["rho"]:+.3f} energy {t["energy"]:+.3f} L_al {t["lal"]:+.3f}  total {t["total"]:+.3f}')

hdr('Q4: points excluded, per region (volume fractions of all points of the region; classes are nested in this order)')
print(f'{"region":15s} {"nonfinite":>9s} {"S<=Smin":>8s} {"cos<=0":>8s} {"k<=0":>6s} {"used":>6s}   cos<=0 share of S>Smin volume')
for rn in ORDER:
    ri = names.index(rn)
    m = reg == ri
    W = w[m].sum()
    fin = np.isfinite(mu1) & np.isfinite(mut) & np.isfinite(c1) & np.isfinite(S1)
    a = m & ~fin
    b = m & fin & ~(S1 > SMIN)
    c = m & fin & (S1 > SMIN) & ~(c1 > 0)
    e = m & fin & (S1 > SMIN) & (c1 > 0) & ~((k1 > 0) & (kr > 0) & (mut > 0) & (mu1 > 0))
    u = m & ok
    sv = m & fin & (S1 > SMIN)
    print(f'{LAB[rn]:15s} {w[a].sum() / W:9.3f} {w[b].sum() / W:8.3f} {w[c].sum() / W:8.3f} {w[e].sum() / W:6.3f} {w[u].sum() / W:6.3f}'
          f'   {w[c].sum() / w[sv].sum():.3f}')

hdr('Q4: bias of the exclusion. B = S>Smin and k>0 and cos>c0 floor; length term needs no cos (L_eff = k^.5|b|/|S| is defined for cos<=0)')
print('terms on mask A (pipeline) and on mask B (every strained point, ln cos floored at c0 = 0.05 / 0.01 for cos < c0):')
for rn in ORDER:
    ri = names.index(rn)
    mA = ok & (reg == ri)
    mB = (reg == ri) & np.isfinite(L1) & np.isfinite(c1) & (S1 > SMIN) & (k1 > 0) & (kr > 0) & (mut > 0) & (L1 > 0)
    out = [f'{LAB[rn]:15s} A: len {mean(ln["length"], mA):+.3f} align {mean(ln["align"], mA):+.3f}']
    out.append(f'B: len {mean(ln["length"], mB):+.3f}')
    for c0 in (0.05, 0.01):
        al = -np.log(np.maximum(c1, c0))
        out.append(f'align(c0={c0}) {mean(al, mB):+.3f}')
    out.append(f'excluded-by-cos share of B volume {1 - w[mA].sum() / w[mB].sum():.3f}')
    print('  '.join(out))

print('\nTail share of the volume mean of the alignment term <-ln cos_193> on mask A:')
print(f'{"region":15s} {"<-lncos>":>8s} {"vol cos<.3":>10s} {"share cos<.3":>12s} {"vol cos<.5":>10s} {"share cos<.5":>12s} {"vol cos>.9":>10s} '
      f'{"share cos>.9":>12s} {"median cos":>10s} {"mean cos":>8s} {"-ln(med)":>8s}')


def wmed(x, ww, q=0.5):
    i = np.argsort(x)
    cw = np.cumsum(ww[i])
    return float(x[i][np.searchsorted(cw, q * cw[-1])])


for rn in ORDER:
    ri = names.index(rn)
    m = ok & (reg == ri)
    al = -np.log(c1)
    tot = (w[m] * al[m]).sum()
    row = f'{LAB[rn]:15s} {tot / w[m].sum():+8.3f} '
    for lo, hi in ((0, 0.3), (0, 0.5), (0.9, 1.01)):
        s = m & (c1 >= lo) & (c1 < hi)
        row += f'{w[s].sum() / w[m].sum():10.3f} {(w[s] * al[s]).sum() / tot:12.3f} '
    row += f'{wmed(c1[m], w[m]):10.3f} {mean(c1, m):8.3f} {-np.log(wmed(c1[m], w[m])):8.3f}'
    print(row)

hdr('Q4: mean of ln against ln of mean, and the gradient-weighted LS regional mu_R = sum(-a:S*)/sum(2|S*|^2) = sum(w mu_eff)/sum(w), w = |S*|^2 dV')
print('(193M mu_eff = rho k^.5 L cos from the interpolated factors; RANS mu_t with its own |S*|_R; RANS has no cos<=0 points)')
print(f'{"region":15s} {"<ln>":>8s} {"ln<>vol":>8s} {"LS|maskA":>9s} {"LS all":>8s} {"LS S>Smin":>9s} {"LS S>Smin,cos>0":>15s}  mu_R193 [Pa s] mu_R_RANS [Pa s]  share of sum(w) with cos<=0 (all)')
fin = np.isfinite(mu1) & np.isfinite(mut) & np.isfinite(S1) & np.isfinite(Sr)
for rn in ORDER:
    ri = names.index(rn)
    mA = ok & (reg == ri)
    dm = lambda x, m: float((dV[m] * x[m]).sum() / dV[m].sum())  # noqa: E731
    lnmean_vol = np.log(dm(mut, mA) / dm(mu1, mA))

    def LS(m):
        m = m & fin
        mr = float((dV[m] * Sr[m] ** 2 * mut[m]).sum() / (dV[m] * Sr[m] ** 2).sum())
        md = float((dV[m] * S1[m] ** 2 * mu1[m]).sum() / (dV[m] * S1[m] ** 2).sum())
        return np.log(mr / md) if (mr > 0 and md > 0) else np.nan, md, mr
    mall = reg == ri
    ms = mall & (S1 > SMIN)
    msc = ms & (c1 > 0)
    la, md_all, mr_all = LS(mall)
    neg = (dV[mall & fin] * S1[mall & fin] ** 2)[(c1[mall & fin] <= 0)].sum() / (dV[mall & fin] * S1[mall & fin] ** 2).sum()
    print(f'{LAB[rn]:15s} {T[rn]["total"]:+8.3f} {lnmean_vol:+8.3f} {LS(mA)[0]:+9.3f} {la:+8.3f} {LS(ms)[0]:+9.3f} {LS(msc)[0]:+15.3f}  '
          f'{md_all:.3e}   {mr_all:.3e}   {neg:.3f}')

# ---------------------------------------------------------------------------------------------- Q3
hdr('Q3: dissipation-length split. 50M 3DNS b006 on its own grid (RANS and 3DNS pointwise, no interpolation), eps = budget eps (resolved)')
rho_d, k_d, L_d, c_d, S_d, me_d, eps = [d[q].astype(np.float64) for q in ('rho_d', 'k_d', 'L_d', 'cos_d', 'S_d', 'meff_d', 'eps')]
with np.errstate(all='ignore'):
    Leps_d = CM * k_d**1.5 / (eps / rho_d)
    Lal_d = me_d / (rho_d * np.sqrt(k_d))
ok50 = (mut > 0) & (me_d > 0) & (c_d > 0) & (S_d > SMIN) & (kr > 0) & (k_d > 0) & (eps > 0) & np.isfinite(Leps_d) & np.isfinite(L_d)
print('RANS L_eps := L_R = mu_t/(rho k^.5) = k^.5/omega (SST, limiter inactive; omega is not exported, so this is an assumption)')
print(f'{"region":15s} {"total50":>8s} | current split: {"dens":>6s} {"energy":>7s} {"length":>7s} {"align":>7s} | L_eps split: {"energy":>7s} {"Leps":>7s} {"struct":>7s} {"comb":>7s} | '
      f'<ln Leff/Leps> <ln Lal/Leps> med(Leff/Leps) med(Lal/Leps)  f_V')
T50 = {}
for rn in ORDER:
    ri = names.index(rn)
    m = ok50 & (reg == ri)
    if m.sum() == 0:
        continue
    with np.errstate(all='ignore'):
        dens = mean(np.log(rr / rho_d), m); en = mean(0.5 * np.log(kr / k_d), m)
        le = mean(np.log(Lr / L_d), m); al = mean(-np.log(c_d), m)
        lee = mean(np.log(Lr / Leps_d), m); st = mean(-np.log(Lal_d / Leps_d), m)
        a1 = mean(np.log(L_d / Leps_d), m); a2 = mean(np.log(Lal_d / Leps_d), m)
        m1 = float(np.exp(wmed(np.log(L_d / Leps_d)[m], w[m]))); m2 = float(np.exp(wmed(np.log(Lal_d / Leps_d)[m], w[m])))
    tot = dens + en + le + al
    T50[rn] = (tot, dens, en, le, al, lee, st)
    print(f'{LAB[rn]:15s} {tot:+8.3f} |               {dens:+6.3f} {en:+7.3f} {le:+7.3f} {al:+7.3f} |              {en:+7.3f} {lee:+7.3f} {st:+7.3f} {lee + st:+7.3f} | '
          f'{a1:+10.3f} {a2:+13.3f} {m1:12.3f} {m2:12.3f}  {w[m].sum() / w[reg == ri].sum():.2f}   chk dens+en+Leps+struct-total = {dens + en + lee + st - tot:+.1e}')

hdr('Q3b: 193M with the proxy eps_193 := eps_50 at the same points (eps of the 193M is not cached). Indicative only.')
with np.errstate(all='ignore'):
    Leps_1 = CM * k1**1.5 / (np.nan_to_num(eps, nan=0) / r1)
    epsok = np.isfinite(eps) & (eps > 0)
print(f'{"region":15s} {"total":>7s} {"density":>7s} {"energy":>7s} {"Leps":>7s} {"struct":>7s}   <ln Leff/Leps> <ln Lal/Leps>')
for rn in ORDER:
    ri = names.index(rn)
    m = ok & epsok & (reg == ri)
    with np.errstate(all='ignore'):
        Lal1 = L1 * c1
        v = [mean(ln['rho'], m), mean(ln['energy'], m), mean(np.log(Lr / Leps_1), m), mean(-np.log(Lal1 / Leps_1), m)]
        a1 = mean(np.log(L1 / Leps_1), m); a2 = mean(np.log(Lal1 / Leps_1), m)
    print(f'{LAB[rn]:15s} {sum(v):+7.3f} {v[0]:+7.3f} {v[1]:+7.3f} {v[2]:+7.3f} {v[3]:+7.3f}   {a1:+12.3f} {a2:+12.3f}   f_V {w[m].sum() / w[reg == ri].sum():.2f}')
