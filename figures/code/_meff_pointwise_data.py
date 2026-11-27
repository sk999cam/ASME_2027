"""Point-wise effective eddy viscosity data (193M b010 against RANS 50M SST), shared by the F14 replacement scripts.

Reads the point file pts.npz (3 141 253 points of the 50M grid with RANS and 193M quantities; written by
planning/maths_audit/audit_extract.py; no 193M read here) and xyz.npz (coordinates, extract_xyz.py).
Mask = "mode C, x/D <= 20" of planning/maths_audit/estimator/estimator.py: finite values, k_193 > 0 and k_RANS > 0
(wall nodes excluded), z/D >= 0.3 on the plate regions, far film restricted to x/D <= 20.
Weights: cell volume V = w/8. All quantities in SI, normalised by rho_inf U_inf D (rho_inf = 0.900 kg/m3,
U_inf = 86.2812 m/s, D = 7.75 mm).
Point-wise estimator: mu_eff = -a:S*/(2|S*|^2) (Wu, Xiao, Paterson 2018), stored in the point file as
rho k^1/2 L cos (193M) and mu_t (RANS).
"""
import numpy as np

SCR = '/tmp/claude-1001/-home-sk999-Documents-ASME-2027/f4c92ceb-6df8-4254-96c3-6552d15d0ecf/scratchpad/'
D, UINF, RHO = 0.00775, 86.2812, 0.88918   # rho_inf = p_exit/(R T_inf)
MUN = RHO * UINF * D
REGIONS = [('near_jet_0_3', 'jet $x/D$ 0–3', 'jet x/D 0-3'), ('wake_3_10', 'near film $x/D$ 3–10', 'near film x/D 3-10'),
           ('wake_10_25', 'far film $x/D$ 10–20', 'far film x/D 10-20')]
QS = (0.10, 0.25, 0.50, 0.75, 0.90)


def wq(x, w, qs=QS):
    """volume-weighted quantiles (midpoint rule on the cumulative weight; negative values kept)."""
    o = np.argsort(x)
    xs, ws = x[o], w[o]
    c = (np.cumsum(ws) - 0.5 * ws) / ws.sum()
    return np.interp(qs, c, xs)


def keep(o, s0, both=False):
    """strain filter: |S*| D/U_inf > s0 for the 193M (both=True: and for RANS)."""
    m = o['SD'] > s0
    return m & (o['SR'] > s0) if both else m


def load():
    d = np.load(SCR + 'pts.npz')
    xyz = np.load(SCR + 'xyz.npz')['xyz']
    names = [str(s) for s in d['regnames']]
    g = lambda n: d[n].astype(np.float64)  # noqa: E731
    reg, V = d['reg'], d['w'].astype(np.float64) / 8.0
    rD, kD, LD, cD, SD = g('c193_rho'), g('c193_k'), g('c193_L'), g('c193_cos'), g('c193_S')
    rR, kR, LR, mut, SR, cR = g('rans_rho'), g('rans_k'), g('rans_L'), g('rans_mut'), g('S_r'), g('rans_cos')
    fin = np.all([np.isfinite(x) for x in (rD, kD, LD, cD, SD, rR, kR, LR, mut, SR)], axis=0)
    valid = fin & (kD > 0) & (kR > 0)
    plate = np.isin(reg, [names.index(r) for r in ('near_jet_0_3', 'wake_3_10', 'wake_10_25', 'upstream_bl')])
    valid &= ~(plate & (xyz[:, 2] < 0.3))
    valid &= ~((reg == names.index('wake_10_25')) & (xyz[:, 0] > 20.0))
    out = {}
    with np.errstate(all='ignore'):
        muD = rD * np.sqrt(kD) * LD * cD
        for key, lab, _ in REGIONS:
            m = valid & (reg == names.index(key))
            sq = lambda a: a[m]  # noqa: E731
            out[key] = dict(
                V=sq(V), muD=sq(muD) / MUN, muR=sq(mut) / MUN,
                kD=sq(kD) / UINF**2, kR=sq(kR) / UINF**2,
                LD=sq(muD) / sq(rD * np.sqrt(kD)) / D, LR=sq(mut) / sq(rR * np.sqrt(kR)) / D,
                cD=sq(cD), cR=sq(cR), SD=sq(SD) * D / UINF, SR=sq(SR) * D / UINF,
                LDfile=sq(LD) / D, LRfile=sq(LR) / D, n=int(m.sum()), nreg=int((reg == names.index(key)).sum()))
    return out
