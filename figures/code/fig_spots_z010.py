#!/usr/bin/env python3
"""fig_spots_z010 (candidate): space-time diagram of the ESTIMATED wall skin friction at mid-pitch, x/D -3 to 20.

C_f,est = 2 mu u(z/D = 0.1) / (0.1 D) / (rho_inf U_inf^2): one-sided gradient from the velocity of the z/D 0.1 plane sampler,
mu the Sutherland viscosity of the local sampler temperature (constants of input_cpu.txt, as asme_n2/n2_extract.py),
rho_inf = p_exit / (R T_inf) (freestream.txt, T_inf = 300 K), U_inf = 86.2812 m/s, D = 7.75 mm.
Line: y/D = +3.00 (sampler j = 300); y/D = -3.00 is identical to float precision (periodic), so no averaging is needed.
Records: every second sampler record of the batch window (the times of the true-wall-C_f set, n2_spots_cf), cached by
_spots_z010_extract.py in /mnt/Elements/DNS/analysis/asme2027_figs/spots_z010/.
(a) 193M 3DNS (C193b, b010) t U/D 150 to 190, (b) 20M LES (C20bW, b006) t U/D 272.8 to 312.8, C_f,est x 10^3, one colour scale
(1st to 99th percentile of the shown data); hole exit x/D -1.0 to 2.35 and the end of the wall data (x/D 14.9) marked.
(c) time-mean of C_f,est over the whole window, time-mean of the true wall C_f (first point above the wall, x/D -3 to 14.88), and
the laminar and turbulent flat-plate correlations of fig_bl_params.py (same Re_D and virtual origin).
Validation (printed, also written to <cache>/validation.txt): C_f,est against the true wall C_f at the same (t, x), x/D -3 to 14.88.
usage: python3 fig_spots_z010.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from cycler import cycler  # noqa: E402

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from functions import set_mpl_params  # noqa: E402
import fig_bl_params as FB  # noqa: E402
import _bl_extract as BX  # noqa: E402
import pandas as pd  # noqa: E402

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_spots_z010.png'
CACHE = '/mnt/Elements/DNS/analysis/asme2027_figs/spots_z010'
N2 = '/mnt/Elements/DNS/analysis/asme2027_figs/n2_spots_cf'
D, U_INF, T_INF, Z = 0.00775, 86.2812, 300.0, 0.1
RUNS = {'c193b': '/mnt/Elements/DNS/film_cooling/C193b_BR15', 'c20bw': '/mnt/Tank/DNS/film_cooling/C20bW_BR15'}
NAME = {'c193b': '193M 3DNS', 'c20bw': '20M LES'}
WIN = {'c193b': (150.0, 190.0), 'c20bw': (272.8, 312.8)}
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
COL = {'c193b': '#000000', 'c20bw': '#D55E00'}
XEND, XHOLE, XWALL = 20.0, (-1.0, 2.35), 14.9
TH_OTSU = 0.00432
GREY = '0.35'


def gas(run):
    """constants of input_cpu.txt (as asme_n2/n2_extract.gas): R, mu_c1, mu_s; p_exit from freestream.txt."""
    for ln in open(os.path.join(run, 'input_cpu.txt')):
        t = ln.split()
        if len(t) == 6:
            try:
                gam, cp, mu_ref, tref, mu_s, prd = map(float, t)
                if 1.0 < gam < 2.0:
                    g = dict(R=cp - cp / gam, mu_c1=mu_ref / tref ** 1.5 * (tref + mu_s), mu_s=mu_s)
                    g['p_exit'] = float(open(os.path.join(run, 'freestream.txt')).read().split()[-1])
                    g['rho_inf'] = g['p_exit'] / (g['R'] * T_INF)
                    return g
            except ValueError:
                pass
    raise SystemExit('no gas line')


def mu_of(T, g):
    return g['mu_c1'] * T.astype(np.float64) ** 1.5 / (T + g['mu_s'])


def load(tag):
    g = gas(RUNS[tag])
    L = np.load(os.path.join(CACHE, f'{tag}_z010_lines.npz'))
    c = np.load(os.path.join(N2, f'{tag}_spacetime.npz'))
    x = L['x'].astype(np.float64)
    u, T, rho = (L[k + '_p'].astype(np.float64) for k in ('u', 'T', 'rho'))
    out = dict(g=g, t=L['t'], x=x, u=u, T=T, rho=rho, mu=mu_of(T, g), muinf=None,
               dm=float(np.abs(L['u_p'] - L['u_m']).max()),
               tc=c['t'], xc=c['x'], cf=c['Cf'].astype(np.float64), cfs=c['Cf_smooth_1p1'].astype(np.float64))
    out['cfe'] = 2.0 * out['mu'] * out['u'] / (Z * D) / (g['rho_inf'] * U_INF ** 2)
    out['cfe_muinf'] = 2.0 * mu_of(np.full(1, T_INF), g)[0] * out['u'] / (Z * D) / (g['rho_inf'] * U_INF ** 2)
    return out


def corr(a, b):
    return float(np.corrcoef(a, b)[0, 1])


def xcorr_lag(a, b, dt, nmax=15):
    """lag (D/U) of the maximum of the cross-correlation of fluctuations; positive: a lags b."""
    a, b = a - a.mean(), b - b.mean()
    n = len(a)
    lags = np.arange(-nmax, nmax + 1)
    r = []
    for k in lags:
        if k >= 0:
            r.append(np.sum(a[k:] * b[:n - k]))
        else:
            r.append(np.sum(a[:n + k] * b[-k:]))
    r = np.array(r) / (n * a.std() * b.std())
    i = int(np.argmax(r))
    if 0 < i < len(r) - 1:
        d = r[i - 1] - 2 * r[i] + r[i + 1]
        sh = 0.5 * (r[i - 1] - r[i + 1]) / d if d != 0 else 0.0
    else:
        sh = 0.0
    return (lags[i] + sh) * dt, float(r[i]), float(r[nmax])


def validate(tag, V, lines):
    """print the validation of C_f,est against the true wall C_f on the overlap."""
    P = lines.append
    xs = V['x'][::2]
    ncx = len(V['xc'])
    est = V['cfe'][:, ::2][:, :ncx]
    tru = V['cf']
    assert np.allclose(xs[:ncx], V['xc'], atol=1e-6)
    sm = V['cfs']
    lam = sm < TH_OTSU
    xc = V['xc']
    P(f'--- {NAME[tag]} ({tag}): {tru.shape[0]} records at the same times, max |dt| sampler-C_f {np.abs(V["t"] - V["tc"]).max():.3f} D/U')
    P(f'    mu: local Sutherland from sampler T (T range {V["T"].min():.2f}..{V["T"].max():.2f} K); mean ratio mu_local/mu_inf(300 K) '
      f'= {np.mean(V["mu"]) / mu_of(np.full(1, T_INF), V["g"])[0]:.4f}; rho_inf {V["g"]["rho_inf"]:.5f} kg/m3; '
      f'max |u(y=+3)-u(y=-3)| {V["dm"]:.2e} m/s')
    P(f'    mean C_f,est(local mu)/C_f,est(mu_inf) over all x, t: {V["cfe"].mean() / V["cfe_muinf"].mean():.4f}')
    # z+ of z/D 0.1 from the true wall C_f
    rho_w, mu_w = V['rho'][:, ::2][:, :ncx], V['mu'][:, ::2][:, :ncx]       # sampler (z/D 0.1) values as wall proxies
    utau = U_INF * np.sqrt(np.maximum(tru, 1e-12) * V['g']['rho_inf'] / (2 * rho_w))
    zp = Z * D * utau * rho_w / mu_w
    bins = [(-3, -1), (2.5, 6), (6, 10), (10, 14.88), (2.5, 14.88), (-3, 14.88)]
    P('    x range        | mean est/true | est/true laminar | est/true turb | corr(est,true) | corr(smoothed) | z+ lam median | z+ turb median | turb fraction')
    for lo, hi in bins:
        m = (xc >= lo) & (xc <= hi)
        e, t_, l = est[:, m], tru[:, m], lam[:, m]
        es = _smooth(est)[:, m]
        z = zp[:, m]
        rl = e[l].mean() / t_[l].mean() if l.any() else np.nan
        rt = e[~l].mean() / t_[~l].mean() if (~l).any() else np.nan
        P(f'    {lo:6.2f}..{hi:6.2f} | {e.mean() / t_.mean():13.3f} | {rl:16.3f} | {rt:13.3f} | {corr(e.ravel(), t_.ravel()):14.3f} | '
          f'{corr(es.ravel(), sm[:, m].ravel()):14.3f} | {np.median(z[l]) if l.any() else np.nan:13.2f} | '
          f'{np.median(z[~l]) if (~l).any() else np.nan:14.2f} | {(~l).mean():.3f}')
    P('    time-mean ratio est/true against x (bins of 1 D, x/D, ratio, corr):')
    for lo in range(-3, 14):
        m = (xc >= lo) & (xc < lo + 1)
        P(f'      {lo:3d}..{lo + 1:3d}  {est[:, m].mean() / tru[:, m].mean():6.3f}  {corr(est[:, m].ravel(), tru[:, m].ravel()):6.3f}')
    P('    lag (cross-correlation of fluctuations, D/U; positive: C_f,est lags the wall C_f), peak r, r at zero lag:')
    dt = float(np.median(np.diff(V['tc'])))
    for xq in (0.0, 4.0, 8.0, 12.0, 14.5):
        j = int(np.argmin(np.abs(xc - xq)))
        lag, rp, r0 = xcorr_lag(est[:, j], tru[:, j], dt)
        P(f'      x/D {xc[j]:6.2f}: lag {lag:6.2f}, r_peak {rp:.3f}, r_0 {r0:.3f}')


def _smooth(a, n=6):
    from scipy.ndimage import uniform_filter1d
    return uniform_filter1d(a, n, axis=0, mode='nearest')


def reference_lines():
    """laminar and turbulent flat-plate C_f x 10^3 against x/D, as fig_bl_params.main (same Re_D and virtual origin)."""
    r50 = pd.read_csv(f'{FB.ST}/{FB.RANS["RANS 50M"]}/post_run1/bl_summary.csv')
    r0 = r50[r50['name'] == 'bl_c_xD-5'].iloc[0]
    nu = r0['ue'] * r0['theta'] / r0['Re_theta']
    ReD = FB.UINF * FB.D / nu
    w50 = BX.load('c50b')
    sel = (np.abs(w50['x'] + 5.0) < 0.15) & (np.abs(w50['y']) < 3.01)
    rex0 = (np.median(w50['ret'][sel]) / 0.664) ** 2
    x0 = -5.0 - rex0 / ReD
    xr = np.linspace(1.0, XEND, 300)
    rex = ReD * (xr - x0)
    return xr, 0.664 / np.sqrt(rex) * 1e3, 0.0576 * rex ** -0.2 * 1e3


def main():
    set_mpl_params(0.5, grid_on=False)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    preview = len(sys.argv) > 1 and sys.argv[1] == 'preview'      # 193M only, separate file
    out = OUT.replace('.png', '_193M_preview.png') if preview else OUT
    tags = [t for t in ('c193b', 'c20bw') if os.path.exists(os.path.join(CACHE, f'{t}_z010_lines.npz'))]
    if preview:
        tags = tags[:1]
    V = {t: load(t) for t in tags}
    lines = []
    for t in tags:
        validate(t, V[t], lines)
    open(os.path.join(CACHE, 'validation.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))

    # colour scale from the data shown
    ix = V[tags[0]]['x'] <= XEND + 1e-9
    shown = []
    for t in tags:
        v = V[t]
        m = (v['t'] >= WIN[t][0]) & (v['t'] <= WIN[t][1])
        shown.append((v['cfe'][m][:, ix] * 1e3).ravel())
    lo, hi = np.percentile(np.concatenate(shown), [1, 99])
    print(f'colour limits (1st, 99th percentile of C_f,est x 1e3): {lo:.3f} {hi:.3f}')

    fig = plt.figure(figsize=(7.2, 8.6))
    nm = len(tags)
    if nm == 1:
        fig.set_size_inches(7.2, 6.0)
    gs = GridSpec(nm + 1, 2, width_ratios=[1, 0.025], height_ratios=[1] * nm + [0.95], hspace=0.12, wspace=0.03,
                  figure=fig, left=0.095, right=0.93, top=0.985, bottom=0.07 if nm == 2 else 0.095)
    axm = []
    im = None
    for i, t in enumerate(tags):
        v = V[t]
        ax = fig.add_subplot(gs[i, 0], sharex=axm[0] if axm else None)
        axm.append(ax)
        m = (v['t'] >= WIN[t][0]) & (v['t'] <= WIN[t][1])
        tt = v['t'][m]
        dtm = np.median(np.diff(tt))
        im = ax.imshow(v['cfe'][m][:, ix] * 1e3, origin='lower', aspect='auto', cmap='inferno', vmin=lo, vmax=hi,
                       extent=[v['x'][0] - 0.01, v['x'][ix][-1] + 0.01, tt[0] - dtm / 2, tt[-1] + dtm / 2],
                       interpolation='nearest')
        ax.set_ylabel(r'$tU_\infty/D$')
        ax.set_xlim(-3, XEND)
        for xh in XHOLE:
            ax.axvline(xh, color='w', ls='--', lw=0.9)
        ax.axvline(XWALL, color='w', ls=':', lw=1.2)
        y0, y1 = ax.get_ylim()
        ax.text(0.5 * sum(XHOLE), y1 - 0.03 * (y1 - y0), 'hole exit', color='w', ha='center', va='top', fontsize=8)
        ax.text(XWALL - 0.25, y1 - 0.03 * (y1 - y0), 'end of wall data', color='w', ha='right', va='top', fontsize=8)
        ax.text(0.012, 0.03, '(%s) %s' % ('ab'[i], NAME[t]), transform=ax.transAxes, color='w', ha='left', va='bottom',
                fontsize=10, bbox=dict(facecolor='k', edgecolor='none', pad=1.5, alpha=0.75))
        plt.setp(ax.get_xticklabels(), visible=False)
    cax = fig.add_subplot(gs[0:nm, 1])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(r'$C_{f,\mathrm{est}}\times10^{3}$')

    ac = fig.add_subplot(gs[nm, 0], sharex=axm[0])
    xr, cl, ct = reference_lines()
    ac.plot(xr, cl, color='0.55', ls='-', lw=1.0)
    ac.plot(xr, ct, color='0.55', ls=':', lw=1.0)
    for t in tags:
        v = V[t]
        ac.plot(v['x'][ix], v['cfe'][:, ix].mean(0) * 1e3, color=COL[t], ls='-')
        ac.plot(v['xc'], v['cf'].mean(0) * 1e3, color=COL[t], ls='--')
    ac.axvline(XWALL, color='k', ls=':', lw=1.0)
    for xh in XHOLE:
        ac.axvline(xh, color='0.5', ls='--', lw=0.8)
    ac.set_xlim(-3, XEND)
    ac.set_xlabel(r'$x/D$')
    ac.set_ylabel(r'$\overline{C_f}\times10^{3}$')
    ac.text(0.012, 0.93, '(c)', transform=ac.transAxes, ha='left', va='top')
    fig.canvas.draw()
    ACF = dict(fontsize=8.5, ha='center', va='center',
               arrowprops=dict(arrowstyle='->', lw=0.7, shrinkA=1, shrinkB=1))
    ylo, yhi = ac.get_ylim()
    ac.set_ylim(ylo, yhi)
    ac._z010 = (xr, cl, ct)
    fig._ac = ac
    annotate_c(ac, V, tags, xr, cl, ct, ACF)
    fig.savefig(out, dpi=300)
    print('saved', out)


def annotate_c(ac, V, tags, xr, cl, ct, ACF):
    """text and arrow labels of panel (c); no legend."""
    ylo, yhi = ac.get_ylim()
    sp = yhi - ylo

    def at(arr_x, arr_y, x):
        return float(np.interp(x, arr_x, arr_y))
    ac.annotate('laminar', xy=(8.0, at(xr, cl, 8.0)), xytext=(8.0, at(xr, cl, 8.0) + 0.14 * sp), color=GREY,
                arrowprops=dict(arrowstyle='->', color=GREY, lw=0.7, shrinkA=1, shrinkB=1), fontsize=8.5, ha='center')
    ac.annotate('turbulent', xy=(6.0, at(xr, ct, 6.0)), xytext=(6.0, at(xr, ct, 6.0) - 0.16 * sp), color=GREY,
                arrowprops=dict(arrowstyle='->', color=GREY, lw=0.7, shrinkA=1, shrinkB=1), fontsize=8.5, ha='center')
    t0 = tags[0]
    v = V[t0]
    ix = v['x'] <= XEND + 1e-9
    kw = dict(arrowprops=dict(arrowstyle='->', lw=0.7, shrinkA=1, shrinkB=1), fontsize=8.5, ha='center')
    ac.annotate('wall, dashed', xy=(11.5, at(v['xc'], v['cf'].mean(0) * 1e3, 11.5)), xytext=(9.0, 0.5 * (ylo + yhi) + 0.25 * sp), **kw)
    ac.annotate(r'est. at $z/D=0.1$, solid', xy=(18.0, at(v['x'][ix], v['cfe'][:, ix].mean(0) * 1e3, 18.0)),
                xytext=(17.0, ylo + 0.12 * sp), **kw)
    for k, t in enumerate(tags):
        ac.text(4.2, ylo + (0.62 - 0.09 * k) * sp, NAME[t], color=COL[t], fontsize=8.5, ha='left', va='center')
    ac.text(XWALL - 0.2, ylo + 0.6 * sp, 'end of wall data', color='k', fontsize=8, ha='right', va='center')
    ac.text(0.5 * sum(XHOLE), ylo + 0.6 * sp, 'hole exit', color='0.4', fontsize=8, ha='center', va='center')


if __name__ == '__main__':
    main()
