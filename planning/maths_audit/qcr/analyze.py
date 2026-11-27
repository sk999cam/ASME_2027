#!/usr/bin/env python3
"""Regional statistics of one cached RANS-50M run (see cache_csv.py).  usage: python3 -B analyze.py TAG [C_CR1]
Writes results_TAG.json next to this script.  Chunked, float64 per chunk; read-only on the cache.
Weights: 'V' = cell volume, 'S2' = V |S*|^2 (S* deviatoric strain from the exported gradients of the same run).
"""
import json
import os
import sys
import numpy as np
import common as cm

tag = sys.argv[1]
CCR = float(sys.argv[2]) if len(sys.argv) > 2 else 0.3
HERE = os.path.dirname(os.path.abspath(__file__))
A1, BS = 0.31, 0.09


def wquant(x, w, qs=(0.1, 0.5, 0.9)):
    o = np.argsort(x)
    c = np.cumsum(w[o])
    c /= c[-1]
    return [float(x[o][min(np.searchsorted(c, q), len(x) - 1)]) for q in qs]


# 1. region membership from X,Y,Z (read once)
X = np.load(f'{cm.C}/{tag}/X.npy', mmap_mode='r')
Y = np.load(f'{cm.C}/{tag}/Y.npy', mmap_mode='r')
Z = np.load(f'{cm.C}/{tag}/Z.npy', mmap_mode='r')
N = len(X)
member = {r: np.zeros(N, bool) for r in cm.REG}
CH = 4_000_000
for i in range(0, N, CH):
    x, y, z = (np.asarray(a[i:i + CH], dtype=np.float64) / cm.D for a in (X, Y, Z))
    for r, f in cm.REG.items():
        member[r][i:i + CH] = f(x, y, z)
res = {'tag': tag, 'c_cr1': CCR, 'N': int(N)}
for reg in cm.REG:
    idx = np.flatnonzero(member[reg])
    d = cm.load(tag, cm.NAMES, idx)
    finite = np.all([np.isfinite(d[n]) for n in cm.NAMES], axis=0)
    ok = finite & (d['k'] > 0) & (d['V'] > 0)                     # wall nodes (k <= 0) and non-finite cells excluded
    out = {'n_cells': int(len(idx)), 'n_used': int(ok.sum()), 'vol_used_D3': float(d['V'][ok].sum() / cm.D ** 3),
           'vol_removed_frac': float(d['V'][~ok].sum() / d['V'].sum())}
    d = {n: v[ok] for n, v in d.items()}
    V = d['V']
    G, S, Ss, W, gn = cm.tensors(d)
    S2 = np.einsum('nij,nij->n', Ss, Ss)
    Sm = np.sqrt(2 * S2)
    wS = V * S2
    # sanity: exported Smag = sqrt(2 S_ij S_ij) (full S) against ours
    Sfull = np.sqrt(2 * np.einsum('nij,nij->n', S, S))
    out['smag_check_max_rel'] = float(np.max(np.abs(Sfull / d['Smag'] - 1)))
    out['divU_over_S_median'] = float(np.median(np.abs(np.trace(G, axis1=1, axis2=2)) / np.maximum(Sfull, 1e-30)))
    mut, rho, k, om = d['mut'], d['rho'], d['k'], d['omega']
    alin = -2 * mut[:, None, None] * Ss
    aq = cm.qcr_part(mut, Ss, W, gn, CCR)           # a = -tau; a_qcr = -tau_NL = +c (O tau - tau O)
    nlin = np.sqrt(np.einsum('nij,nij->n', alin, alin))
    nq = np.sqrt(np.einsum('nij,nij->n', aq, aq))
    nS = np.sqrt(S2)
    cs = np.einsum('nij,nij->n', aq, Ss) / np.maximum(nq * nS, 1e-300)
    m = nq > 0
    out['qcr_contraction_max'] = float(np.max(np.abs(cs[m])))
    out['qcr_trace_max_rel'] = float(np.max(np.abs(np.trace(aq, axis1=1, axis2=2))[m] / nq[m]))
    out['qcr_antisym_max_rel'] = float(np.max(np.abs(aq - np.swapaxes(aq, 1, 2)).max(axis=(1, 2))[m] / nq[m]))
    ratio = nq / np.maximum(nlin, 1e-300)
    cos = nlin / np.sqrt(nlin ** 2 + nq ** 2)
    for wn, w in (('V', V), ('S2', wS)):
        sw = w.sum()
        q = {}
        q['ratio_agg'] = float(np.sqrt((w * nq ** 2).sum() / (w * nlin ** 2).sum()))
        q['ratio_q10_50_90'] = wquant(ratio, w)
        q['cos_mean'] = float((w * cos).sum() / sw)
        q['cos_q10_50_90'] = wquant(cos, w)
        # a_perp components in plate axes (x streamwise, y lateral, z wall-normal), normalised by <rho k>_w
        rk = (w * rho * k).sum() / sw
        A = {}
        for i, a in enumerate('xyz'):
            for j, b in enumerate('xyz'):
                if j >= i:
                    A[a + b] = float((w[:, ] * aq[:, i, j]).sum() / sw / rk)
        q['aperp_over_rhok_mean'] = A
        q['aperp_yy-zz'] = A['yy'] - A['zz']
        q['aperp_xx-yy'] = A['xx'] - A['yy']
        q['aperp_xx-zz'] = A['xx'] - A['zz']
        # rms of the same components, for the sign-cancellation check
        q['aperp_rms_over_rhok'] = {a + b: float(np.sqrt((w * aq[:, i, j] ** 2).sum() / sw) / rk)
                                    for i, a in enumerate('xyz') for j, b in enumerate('xyz') if j >= i}
        # linear-part components for reference
        q['alin_over_rhok_mean'] = {a + b: float((w * alin[:, i, j]).sum() / sw / rk)
                                    for i, a in enumerate('xyz') for j, b in enumerate('xyz') if j >= i}
        out['qcr_stress_' + wn] = q
    # task 3: mean-field quantities of the run itself (own gradients)
    sV, sS = V.sum(), wS.sum()
    out['muR'] = float((wS * mut).sum() / sS)                              # LS: sum V mu_t |S*|^2 / sum V |S*|^2
    out['mut_vol_mean'] = float((V * mut).sum() / sV)
    out['k_vol_mean'] = float((V * k).sum() / sV)
    out['rhok_vol_mean'] = float((V * rho * k).sum() / sV)
    out['k_S2_mean'] = float((wS * k).sum() / sS)
    wx = d['dWdy'] - d['dVdz']
    out['omega_x_rms'] = float(np.sqrt((V * wx ** 2).sum() / sV))
    out['omega_x_rms_over_U_D'] = out['omega_x_rms'] * cm.D / cm.U_INF
    out['v_rms'] = float(np.sqrt((V * d['v'] ** 2).sum() / sV))
    out['w_rms'] = float(np.sqrt((V * d['w'] ** 2).sum() / sV))
    out['v_mean'] = float((V * d['v']).sum() / sV)
    out['dTdy_sq_int'] = float((V * d['dTdy'] ** 2).sum() / cm.D ** 3)     # K^2/m^2 * m^3 / D^3 -> K^2 /m^2 * (volume in D^3)
    out['dTdy_rms'] = float(np.sqrt((V * d['dTdy'] ** 2).sum() / sV))
    out['dTdz_sq_int'] = float((V * d['dTdz'] ** 2).sum() / cm.D ** 3)
    out['T_vol_mean'] = float((V * d['T']).sum() / sV)
    out['S_rms'] = float(np.sqrt((V * S2).sum() / sV))
    # LS viscosity with the full stress of this run: -a:S*/(2|S*|^2); QCR part contributes nothing
    a_full = alin + aq
    out['muR_full_stress'] = float((V * (-np.einsum('nij,nij->n', a_full, Ss))).sum() / (2 * sS))
    out['muR_qcr_part_contrib'] = float((V * (-np.einsum('nij,nij->n', aq, Ss))).sum() / (2 * sS))
    # task 4: SST limiter / dissipation length
    wd, mu = d['wd'], d['mu']
    nu = mu / rho
    arg2 = np.maximum(2 * np.sqrt(k) / (BS * om * wd), 500 * nu / (wd ** 2 * om))
    F2 = np.tanh(arg2 ** 2)
    Lq = np.sqrt(k) / om                                       # C_mu k^1.5/eps with eps = beta* k omega, C_mu = beta*
    Lmu = mut / (rho * np.sqrt(k))
    r = Lmu / Lq                                               # = omega mu_t/(rho k) = omega/omega_tilde
    r_pred = np.minimum(1.0, A1 * om / np.maximum(Sm * F2, 1e-300))
    act = r < 0.999
    out['lim'] = {
        'ratio_mu_t_over_rho_sqrtk_L_eps_V_mean': float((V * r).sum() / sV),
        'ratio_S2_mean': float((wS * r).sum() / sS),
        'ratio_V_q10_50_90': wquant(r, V),
        'ratio_S2_q10_50_90': wquant(r, wS),
        'ratio_LS_regional': float((wS * rho * np.sqrt(k) * Lq * r).sum() / (wS * rho * np.sqrt(k) * Lq).sum()),
        'active_vol_frac': float(V[act].sum() / sV),
        'active_S2_frac': float(wS[act].sum() / sS),
        'L_eps_over_D_V_median': wquant(Lq / cm.D, V, (0.5,))[0],
        'L_mu_over_D_V_median': wquant(Lmu / cm.D, V, (0.5,))[0],
        'L_eps_S2_mean_over_D': float((wS * Lq).sum() / sS / cm.D),
        'L_mu_S2_mean_over_D': float((wS * Lmu).sum() / sS / cm.D),
        'F2_V_median': wquant(F2, V, (0.5,))[0],
        'F2_S2_mean': float((wS * F2).sum() / sS),
        'pred_vs_actual_ratio_max_abs': float(np.max(np.abs(r_pred - r))),
        'pred_vs_actual_ratio_p99_abs': float(np.quantile(np.abs(r_pred - r), 0.99)),
        'pred_active_agree_frac': float(np.mean((r_pred < 0.999) == act)),
        'r_max': float(r.max()), 'r_min': float(r.min()),
    }
    # sanity ranges
    out['sanity'] = {'k_min': float(k.min()), 'k_max': float(k.max()), 'mut_over_mu_max': float((mut / mu).max()),
                     'rho_min': float(rho.min()), 'rho_max': float(rho.max()), 'T_min': float(d['T'].min()), 'T_max': float(d['T'].max())}
    res[reg] = out
    print(tag, reg, "done", out['n_used'], flush=True)
json.dump(res, open(os.path.join(HERE, f'results_{tag}.json'), 'w'), indent=1)
