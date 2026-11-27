"""Near-wall / far-wall split of the eddy-viscosity difference for fig_meff_maps.py (helper module).

Logic of split_compute.py (review session), without file output. Reads the point file through _meff_pointwise_data.SCR.
Region value (route A): mu_R = sum V(-a:S*) / sum 2V|S*|^2 on each side with its own |S*|; ln(mu_R,RANS/mu_R,193M) is split
exactly (Shapley, all orders) over rho, k^1/2, L_eff, cos(theta), strain weights.
Sub-regions: near film x/D 3-10 and far film x/D 10-20; near wall 0 < z < z_x(x), far wall z_x <= z <= 2 D; wall nodes excluded.
"""
import itertools

import numpy as np

import _meff_pointwise_data as M

D, UINF, RHO, MUN = M.D, M.UINF, M.RHO, M.MUN
TERMS = ['rho', 'k12', 'Leff', 'cos', 'S']
_C = {}


def _load():
    if _C:
        return _C
    SCR = M.SCR
    d = np.load(SCR + 'pts.npz')
    X = np.load(SCR + 'xyz.npz')['xyz']
    names = [str(s) for s in d['regnames']]
    g = lambda n: d[n].astype(np.float64)  # noqa: E731
    reg, V = d['reg'], d['w'].astype(np.float64) / 8.0
    q = dict(rD=g('c193_rho'), kD=g('c193_k'), LD=g('c193_L'), cD=g('c193_cos'), SD=g('c193_S'),
             rR=g('rans_rho'), kR=g('rans_k'), LR=g('rans_L'), mut=g('rans_mut'), SR=g('S_r'))
    fin = np.all([np.isfinite(v) for v in q.values()], axis=0)
    base = fin & (q['kD'] > 0) & (q['kR'] > 0)
    with np.errstate(all='ignore'):
        LRd = q['mut'] / (q['rR'] * np.sqrt(q['kR']))                       # RANS L_R = mu_t / (rho k^1/2)
        muD = q['rD'] * np.sqrt(q['kD']) * q['LD'] * q['cD']                 # 193M mu_eff = rho k^1/2 L_eff cos
    _C.update(q=q, X=X, names=names, reg=reg, V=V, base=base, LRd=LRd, muD=muD)
    return _C


def _sums(c, m):
    q, V = c['q'], c['V'][m]
    A = dict(rho=(q['rD'][m], q['rR'][m]), k=(np.sqrt(q['kD'][m]), np.sqrt(q['kR'][m])), L=(q['LD'][m], q['LR'][m]),
             c=(q['cD'][m], np.ones(m.sum())), S=(q['SD'][m], q['SR'][m]))
    keys = ['rho', 'k', 'L', 'c', 'S']
    fv = {}
    for s in itertools.product([False, True], repeat=5):
        xx = {k: A[k][int(s[i])] for i, k in enumerate(keys)}
        den = V * xx['S'] ** 2
        r = (den * xx['rho'] * xx['k'] * xx['L'] * xx['c']).sum() / den.sum()
        fv[s] = np.log(r) if r > 0 else np.nan
    return fv


def _shapley(fv):
    phi = np.zeros(5)
    ords = list(itertools.permutations(range(5)))
    for o in ords:
        st = [False] * 5
        prev = fv[tuple(st)]
        for i in o:
            st[i] = True
            cur = fv[tuple(st)]
            phi[i] += (cur - prev) / len(ords)
            prev = cur
    return phi


def _wmed(a, w):
    return M.wq(a, w, (0.5,))[0]


def compute(XT, ZT, ztop=2.0, layers=None):
    """Returns {(region, layer): dict} for region in near/far film.
    Default layers: near wall 0 <= z < z_x(x), far wall z_x <= z <= ztop, z_x interpolated linearly from the table (XT, ZT).
    layers: {name: (lo, hi)}, each bound None (wall / ztop) or a table (XT, ZT); lower bound inclusive, upper exclusive."""
    c = _load()
    q, X, V, reg, names, base = c['q'], c['X'], c['V'], c['reg'], c['names'], c['base']
    x, z = X[:, 0], X[:, 2]
    if layers is None:
        layers = {'near wall': (None, (XT, ZT)), 'far wall': ((XT, ZT), None)}
    regs = {'near film': (x >= 3) & (x < 10) & (reg == names.index('wake_3_10')),
            'far film': (x >= 10) & (x <= 20) & (reg == names.index('wake_10_25'))}
    out = {}
    for rk, rm in regs.items():
        b = base & rm & (z <= ztop + 1e-9)
        Vt = V[b].sum()
        WD = (V * q['SD'] ** 2)[b].sum()
        for lk, (lo, hi) in layers.items():
            lm = np.ones(len(z), bool)
            if lo is not None:
                lm &= z >= np.interp(x, *lo)
            if hi is not None:
                lm &= z < np.interp(x, *hi)
            m = b & lm
            fv = _sums(c, m)
            p5 = _shapley(fv)
            tot = fv[(True,) * 5] - fv[(False,) * 5]
            v = V[m]
            muD = q['rD'][m] * np.sqrt(q['kD'][m]) * q['LD'][m] * q['cD'][m]
            muR = q['mut'][m]
            r = dict(tot=tot, p5=p5, res=p5.sum() - tot,
                     plain=np.log((v * muR).sum() / (v * muD).sum()),
                     vshare=v.sum() / Vt, sshare=(V * q['SD'] ** 2)[m].sum() / WD,
                     fneg=v[muD <= 0].sum() / v.sum(), n=int(m.sum()))
            with np.errstate(all='ignore'):
                for tag, sel in (('all', np.ones(m.sum(), bool)), ('S05', q['SD'][m] * D / UINF > 0.05)):
                    r['med_' + tag] = _wmed(muR[sel], v[sel]) / _wmed(muD[sel], v[sel])
                    r['keep_' + tag] = v[sel].sum() / v.sum()
            out[(rk, lk)] = r
    return out


def _shapley4(fv):
    """exact Shapley split over four players (rho, k^1/2, L_eff, cos); fv maps a 4-tuple of booleans to ln of the region sum."""
    phi = np.zeros(4)
    ords = list(itertools.permutations(range(4)))
    for o in ords:
        st = [False] * 4
        prev = fv[tuple(st)]
        for i in o:
            st[i] = True
            cur = fv[tuple(st)]
            phi[i] += (cur - prev) / len(ords)
            prev = cur
    return phi


def split_mask(m):
    """Plain volume statistics over the points selected by the boolean mask m: total = ln(sum V mu_t,RANS / sum V mu_eff,193M), exact Shapley
    split into rho, k^1/2, L_eff, cos (four players), the volume-weighted means of each factor for both rungs, and the share with mu_eff <= 0."""
    c = _load()
    q, V, LRd, muD = c['q'], c['V'], c['LRd'], c['muD']
    v = V[m]
    A = dict(rho=(q['rD'][m], q['rR'][m]), k=(np.sqrt(q['kD'][m]), np.sqrt(q['kR'][m])), L=(q['LD'][m], LRd[m]),
             c=(q['cD'][m], np.ones(m.sum())))
    keys = ['rho', 'k', 'L', 'c']
    fv = {}
    for s in itertools.product([False, True], repeat=4):
        prod = np.ones(m.sum())
        for i, kk in enumerate(keys):
            prod = prod * A[kk][int(s[i])]
        r = (v * prod).sum()
        fv[s] = np.log(r) if r > 0 else np.nan
    tot = np.log((v * q['mut'][m]).sum() / (v * muD[m]).sum())
    p4 = _shapley4(fv)
    vm = {kk: ((v * A[kk][1]).sum() / v.sum(), (v * A[kk][0]).sum() / v.sum()) for kk in keys}       # (RANS, 193M) volume means
    return dict(tot=tot, p4=p4, res=p4.sum() - tot, chk=(fv[(True,) * 4] - fv[(False,) * 4]) - tot, vmean=vm,
                fneg=v[muD[m] <= 0].sum() / v.sum(), n=int(m.sum()))


def compute_curve(zx_tab, centres, width=1.0, s0=0.05, ztop=2.0, layer='near'):
    """split_mask in slabs of x (|x - centre| < width/2) over the whole film x/D 3-20, for the layer below z_x(x) (layer 'near': wall nodes
    excluded, z < z_x) or above it ('above': z >= z_x up to the strain-mask edge); z_x(x) from the table zx_tab = (x, z_x), interpolated
    linearly at each point. Mask: point-wise 193M |S*| D/U_inf > s0. Returns a list of dicts (split_mask result plus x, kept share of the layer
    volume, mean z of the kept points)."""
    c = _load()
    X, V, reg, names, base, q = c['X'], c['V'], c['reg'], c['names'], c['base'], c['q']
    x, z = X[:, 0], X[:, 2]
    film = ((x >= 3) & (x < 10) & (reg == names.index('wake_3_10'))) | ((x >= 10) & (x <= 20) & (reg == names.index('wake_10_25')))
    b = base & film & (z <= ztop + 1e-9)
    zx = np.interp(x, *zx_tab)
    lay = (z < zx) if layer == 'near' else (z >= zx)
    sd = q['SD'] * D / UINF
    out = []
    for xc_ in np.atleast_1d(centres):
        slab = b & lay & (np.abs(x - xc_) < width / 2)
        m = slab & (sd > s0)
        r = split_mask(m)
        r.update(x=float(xc_), kept=V[m].sum() / V[slab].sum(), zmean=float((V[m] * z[m]).sum() / V[m].sum()))
        out.append(r)
    return out


def compute_mask(s0=0.05, ztop=2.0, zsplit=0.6):
    """Strain as a mask only: cells with the point-wise 193M |S*| D/U_inf > s0 are kept, everything inside is plain volume-weighted.
    Regions near film x/D 3-10 and far film 10-20, z <= ztop, wall nodes excluded (k_193 = 0 or k_RANS <= 0 or non-finite).
    Primary: total = ln(sum V mu_t,RANS / sum V mu_eff,193M), exact Shapley split over rho, k^1/2, L_eff, cos (RANS: cos = 1,
    L_R = mu_t / (rho k^1/2)). Secondary: volume mean of the point-wise ln ratios (cells with mu_eff,193M > 0) and ratios of
    volume-weighted medians. Returns {(region, band): dict}, band in all / z<0.6 / z>=0.6."""
    c = _load()
    q, X, V, reg, names, base = c['q'], c['X'], c['V'], c['reg'], c['names'], c['base']
    x, z = X[:, 0], X[:, 2]
    with np.errstate(all='ignore'):
        LRd = q['mut'] / (q['rR'] * np.sqrt(q['kR']))
        muD = q['rD'] * np.sqrt(q['kD']) * q['LD'] * q['cD']
    regs = {'near film': (x >= 3) & (x < 10) & (reg == names.index('wake_3_10')),
            'far film': (x >= 10) & (x <= 20) & (reg == names.index('wake_10_25'))}
    regs['whole film'] = regs['near film'] | regs['far film']          # x/D 3-20 merged
    sd = q['SD'] * D / UINF
    out = {}
    for rk, rm in regs.items():
        b = base & rm & (z <= ztop + 1e-9)
        keep = b & (sd > s0)
        for bk, bm in (('all', np.ones(len(z), bool)), ('z<0.6', z < zsplit), ('z>=0.6', z >= zsplit)):
            m = keep & bm
            v = V[m]
            core = split_mask(m)
            tot, p4, chk = core['tot'], core['p4'], core['chk']
            pos = muD[m] > 0
            vp = v[pos]
            lnr = dict(rho=np.log(q['rR'][m][pos] / q['rD'][m][pos]), k=np.log(np.sqrt(q['kR'][m][pos] / q['kD'][m][pos])),
                       L=np.log(LRd[m][pos] / q['LD'][m][pos]), c=-np.log(q['cD'][m][pos]))
            lnmean = {k_: (vp * a).sum() / vp.sum() for k_, a in lnr.items()}
            lnmean['tot'] = sum(lnmean.values())
            med = {}
            for key, (aR, aD) in dict(rho=(q['rR'][m], q['rD'][m]), k=(np.sqrt(q['kR'][m]), np.sqrt(q['kD'][m])), L=(LRd[m], q['LD'][m]),
                                      c=(np.ones(m.sum()), q['cD'][m]), mu=(q['mut'][m], muD[m])).items():
                med[key] = _wmed(aR, v) / _wmed(aD, v)
            xb = np.floor((x[m] - 3.0) / 0.25).astype(int)
            top = np.array([z[m][xb == i].max() for i in np.unique(xb)]) if m.any() else np.array([np.nan])
            Vb = V[b & bm].sum()
            out[(rk, bk)] = dict(tot=tot, p4=p4, res=p4.sum() - tot, chk=chk, kept=V[m].sum() / Vb, keptall=V[keep].sum() / V[b].sum(),
                                 fneg=v[~pos].sum() / v.sum(), lnmean=lnmean, med=med, top=float(np.mean(top)), n=int(m.sum()),
                                 plain=tot, muRm=float((v * q['mut'][m]).sum() / v.sum() / M.MUN),
                                 muDm=float((v * muD[m]).sum() / v.sum() / M.MUN))
    return out
