"""Shared definitions for the QCR analysis (read-only on inputs; arrays cached by cache_csv.py)."""
import numpy as np

C = '/mnt/Elements/DNS/analysis/asme2027_figs/qcr'
D = 0.00775                      # hole diameter [m]: dns/python/make_tke_inputs.py line 43; analysis/meff_c50b_c100b/errorbars.py
U_INF = 86.2812
ZCUT = 0.3
# Region bounds = dns/python/make_tke_inputs.py REGIONS (via analysis/regions_pp/build_regions.py), applied to cell centres
# with z/D >= 0.3 (estimator mode C).  jet: 0<=x<=3, |y|<=1.5, z<=2 ; film 3-10: 3<x<=10, z<=2 ; film 10-20: 10<x<=20, z<=2.
# (no |y| bound on the film regions, as in make_tke_inputs; the 50M domain is one pitch).
REG = {
    'jet 0-3': lambda x, y, z: (z >= ZCUT) & (z <= 2) & (x >= 0) & (x <= 3) & (np.abs(y) <= 1.5),
    'film 3-10': lambda x, y, z: (z >= ZCUT) & (z <= 2) & (x > 3) & (x <= 10),
    'film 10-20': lambda x, y, z: (z >= ZCUT) & (z <= 2) & (x > 10) & (x <= 20),
    # indicative centre-plane band (not an estimator region): film 3-20, |y|/D <= 0.25
    'centre band 3-20': lambda x, y, z: (z >= ZCUT) & (z <= 2) & (x > 3) & (x <= 20) & (np.abs(y) <= 0.25),
}
NAMES = ['u', 'v', 'w', 'dUdx', 'dUdy', 'dUdz', 'dVdx', 'dVdy', 'dVdz', 'dWdx', 'dWdy', 'dWdz', 'Smag', 'T',
         'dTdx', 'dTdy', 'dTdz', 'rho', 'mu', 'k', 'omega', 'mut', 'wd', 'gamma', 'V', 'X', 'Y', 'Z']


def load(tag, names, idx=None):
    out = {}
    for n in names:
        a = np.load(f'{C}/{tag}/{n}.npy', mmap_mode='r')
        out[n] = np.asarray(a if idx is None else a[idx], dtype=np.float64)
    return out


def tensors(d):
    """Return G (n,3,3) = dU_i/dx_j, S* (deviatoric strain), W (rotation), |grad v|."""
    n = len(d['u'])
    G = np.empty((n, 3, 3))
    for i, c in enumerate('UVW'):
        for j, x in enumerate('xyz'):
            G[:, i, j] = d[f'd{c}d{x}']
    S = 0.5 * (G + np.swapaxes(G, 1, 2))
    W = 0.5 * (G - np.swapaxes(G, 1, 2))
    tr = np.trace(S, axis1=1, axis2=2)
    Ss = S.copy()
    for i in range(3):
        Ss[:, i, i] -= tr / 3
    gn = np.sqrt(np.einsum('nij,nij->n', G, G))
    return G, S, Ss, W, gn


def qcr_part(mut, Ss, W, gn, c_cr1=0.3):
    """a_qcr = -tau_NL (a = -tau, tau = Boussinesq stress incl. sign of STAR's T_RANS), standard QCR2000:
    tau_NL = -c_cr1 (O_ik tau_jk + O_jk tau_ik), tau = 2 mu_t S*, O = 2 W/|grad v|.
    O_ik tau_jk + O_jk tau_ik = (O tau - tau O)_ij  (tau symmetric, O antisymmetric)."""
    O = 2 * W / gn[:, None, None]
    tau = 2 * mut[:, None, None] * Ss
    comm = np.einsum('nik,nkj->nij', O, tau) - np.einsum('nik,nkj->nij', tau, O)
    return c_cr1 * comm            # a_qcr = + c_cr1 (O tau - tau O)
