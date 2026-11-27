#!/usr/bin/env python3
"""Q5: does the QCR correction have zero contraction with S*?

QCR2000 (Spalart 2000), tau = -rho u''u'' (so a = -tau):
    tau_QCR_ij = tau_lin_ij - c_cr1 (O_ik tau_jk + O_jk tau_ik),  O_ik = 2 Omega_ik / sqrt(dU_m/dx_n dU_m/dx_n),
    Omega_ik = 1/2 (dU_i/dx_k - dU_k/dx_i),  tau_lin = 2 mu_t S* - 2/3 rho k I  (sign of tau: tau = -2 mu_t S* in the
    task's wording means a = +2 mu_t S*? No: a = -tau_... see below; the identity does not depend on the sign).
Index form: O_ik tau_jk + O_jk tau_ik = (O tau^T + tau O^T)... with tau symmetric = O tau - tau O = [O, tau].
With tau = c S* + d I: [O, tau] = c [O, S*]  (I commutes), and  [O,S*]:S* = tr(O S* S*) - tr(S* O S*) = 0 by cyclicity.
This script checks that numerically (random tensors), builds a_perp = a + 2 mu_eff S* for a QCR stress, and gives
cos(theta) of the QCR stress in simple shear.
usage: python3 -B q5_qcr.py
"""
import numpy as np

rng = np.random.default_rng(7)
N = 100000


def qcr_stress(G, mu_t, rho_k, c=0.307, ccr=True):
    """tau_ij = -rho u''_i u''_j (Reynolds stress with sign tau) of QCR2000 for gradient G_ij = dU_i/dx_j."""
    S = 0.5 * (G + np.swapaxes(G, -1, -2))
    Om = 0.5 * (G - np.swapaxes(G, -1, -2))
    tr = np.trace(S, axis1=-2, axis2=-1)
    Ss = S.copy()
    for i in range(3):
        Ss[:, i, i] -= tr / 3
    Gn = np.sqrt(np.einsum('nij,nij->n', G, G))
    O = 2 * Om / Gn[:, None, None]
    tau_lin = 2 * mu_t[:, None, None] * Ss - 2.0 / 3.0 * rho_k[:, None, None] * np.eye(3)
    # explicit index form of the task: O_ik tau_jk + O_jk tau_ik
    corr = np.einsum('nik,njk->nij', O, tau_lin) + np.einsum('njk,nik->nij', O, tau_lin)
    return Ss, O, tau_lin, (tau_lin - c * corr if ccr else tau_lin), corr


G = rng.normal(0, 1e3, (N, 3, 3))
mu_t = rng.uniform(1e-3, 1e-1, N)
rho_k = rng.uniform(1, 50, N)
Ss, O, tau_lin, tau, corr = qcr_stress(G, mu_t, rho_k)
S2 = np.einsum('nij,nij->n', Ss, Ss)
cs = np.einsum('nij,nij->n', corr, Ss)
cn = np.sqrt(np.einsum('nij,nij->n', corr, corr))
print('QCR correction c:S* / (|c||S*|): max |.| = %.2e (n = %d random tensors)' % (np.abs(cs / (cn * np.sqrt(S2))).max(), N))
print('correction symmetric: max |c - c^T|/|c| = %.2e ; trace/|c| max = %.2e' % (
    (np.abs(corr - np.swapaxes(corr, -1, -2)).max(axis=(1, 2)) / cn).max(), (np.abs(np.trace(corr, axis1=1, axis2=2)) / cn).max()))
# mu_eff of the QCR stress: a = -tau (deviatoric part); mu_eff = -a:S*/(2|S*|^2) = tau:S*/(2|S*|^2)
dev = lambda T: T - (np.trace(T, axis1=1, axis2=2) / 3)[:, None, None] * np.eye(3)  # noqa: E731
a = -dev(tau)
mu_eff = -np.einsum('nij,nij->n', a, Ss) / (2 * S2)
print('mu_eff(QCR)/mu_t - 1: max |.| = %.2e' % np.abs(mu_eff / mu_t - 1).max())
aperp = a + 2 * mu_eff[:, None, None] * Ss
an = np.sqrt(np.einsum('nij,nij->n', a, a))
ap = np.sqrt(np.einsum('nij,nij->n', aperp, aperp))
print('a_perp = a + 2 mu_eff S*: |a_perp|/|a| median %.3f ; a_perp traceless: max |tr|/|a_perp| = %.1e' % (
    np.median(ap / an), (np.abs(np.trace(aperp, axis1=1, axis2=2)) / ap).max()))
# simple shear U = gamma z: G_13 = gamma (u_x gradient in z), others 0
gam = 1e4
G1 = np.zeros((1, 3, 3)); G1[0, 0, 2] = gam
for c in (0.3, 0.307):
    Ss1, O1, tl, t1, c1 = qcr_stress(G1, np.array([0.01]), np.array([10.0]), c=c)
    a1 = -dev(t1); S21 = np.einsum('nij,nij->n', Ss1, Ss1)
    me = -np.einsum('nij,nij->n', a1, Ss1) / (2 * S21)
    ap1 = a1 + 2 * me[:, None, None] * Ss1
    r = np.sqrt(np.einsum('nij,nij->n', ap1, ap1) / np.einsum('nij,nij->n', a1 - ap1, a1 - ap1))
    print('simple shear, c_cr1 = %.3f: mu_eff/mu_t = %.6f ; |a_perp|/|a_par| = %.4f (= 2 c_cr1) ; cos(theta) = %.4f ; ln cos = %.4f' % (
        c, me[0] / 0.01, r[0], 1 / np.sqrt(1 + r[0]**2), np.log(1 / np.sqrt(1 + r[0]**2))))
    print('   a_perp (units of mu_t gamma):\n', np.round(ap1[0] / (0.01 * gam), 4))
