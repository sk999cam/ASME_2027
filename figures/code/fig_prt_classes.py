"""fig_prt_classes: classic turbulent Prandtl number per flow class (least squares), 193M reference.

Shows Pr_t,R = [sum(-rho u''w'' dU/dz) / sum(dU/dz)^2] / [sum(-rho w''T'' dT/dz) / sum(dT/dz)^2] (Favre means,
wall-normal classic form of Kohli and Bogard 2005) over classes of points of the 193M 3DNS (C193b): stat3 batch b010
(10 FT cumulative mean, filled markers, bars = 95 % block-bootstrap interval) and b006 (6 FT, hollow markers).
No threshold: a point enters with the weight of its gradient squared (ratio of sums, not mean of ratios), so a
near-zero gradient cannot produce a near-zero denominator. The tensor form Pr_C = mu_C/alpha_C (c_p mu_eff/alpha_eff,
same estimator) is drawn as small grey markers (b010) for the methods-appendix comparison.
Classes: (a) wall-distance bands z/D, per region (jet x/D 0-3, near film 3-10, far film 10-20 -> 20 here) and the
cross-plane x/D = 10; (b) coolant-concentration bands eta (z/D >= 0.3; the wall layer z/D < 0.3 is a class of its
own); (c) fit quality of a scalar diffusivity per z/D band: R2_flux = 1 - sum|rho q + alpha_C grad T|^2 / sum|rho q|^2
(filled) and R2_stress = 1 - sum|a + 2 mu_C S*|^2 / sum|a|^2 (hollow), a = rho R minus its trace part.
Lines: RANS Pr_t = 0.9; Kohli and Bogard (2005, J. Heat Transfer 127:513), DR 1.05, M = 0.4: about 2 at low
freestream turbulence (TL 0.6 %), 0.9 to 1.2 at high freestream turbulence (TL 20 %).
Bootstrap: points are grouped in 1 D wide x bins (centre plane) or 0.5 D wide y bins (cross-plane), the bins are
resampled with replacement (1000 draws), the 2.5 and 97.5 percentiles are drawn. Values outside the y range are not
drawn; the tables printed by the script and written to the cache folder give every value.
Data: /mnt/Elements/DNS/film_cooling/C193b_BR15, stat3_b006 and stat3_b010, through the plane caches of
meff_planes_common.py (tensor form) and prt_corr_common.py (components; whole centre plane x/D 0-20, cross-plane
x/D = 10). Cache: /mnt/Elements/DNS/analysis/asme2027_figs/prt_classes/. Complete, not interim, when both b006 and
b010 caches exist.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.lines import Line2D

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, '/home/sk999/Documents/ASME_2027/figures/code')
from functions import set_mpl_params
import meff_planes_common as mc
import prt_corr_common as cc

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_prt_classes.pdf'
TABLE = os.path.join(cc.CORR_DIR, 'prt_classes_tables.txt')
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
D = mc.D
NBOOT = 1000
REGIONS = [('jet', 0.0, 3.0), ('near film', 3.0, 10.0), ('far film', 10.0, 20.0)]
ZB = [0.0, 0.05, 0.15, 0.3, 0.6, 3.01]
ZLAB = ['0-0.05', '0.05-0.15', '0.15-0.3', '0.3-0.6', '>0.6']
EB = [('core', 0.4, 9.0), ('mixing layer', 0.1, 0.4), ('edge', 0.02, 0.1), ('freestream', -9.0, 0.02)]
KEYS = ['nm', 'dm', 'na', 'da', 'a2', 'cnm', 'cdm', 'cna', 'cda', 'lnm', 'ldm', 'lna', 'lda', 'F2', 'g2', 'Fg']
LINES = []


def log(s=''):
    print(s)
    LINES.append(s)


def load_corr(tag, batch, shape):
    z = dict(np.load(os.path.join(cc.CORR_DIR, f'corr_{tag}_b{batch:03d}.npz')))
    out = {}
    for q in cc.FIELDS:
        a = np.full(int(np.prod(shape)), np.nan, np.float32)
        a[z['idx']] = z[q]
        out[q] = a.astype(float)
    return out


def prep(R, c, p):
    """Flat per-point arrays: coordinates, bins, class variables, sums terms (zero where undefined)."""
    S, G = R[f'{p}_S'], R[f'{p}_G']
    me, ae = R[f'{p}_meff'].ravel(), R[f'{p}_aeff'].ravel()
    S, G = S.ravel(), G.ravel()
    A, B = R[f'{p}_A'].ravel() / D, R[f'{p}_B'].ravel() / D
    rho = R[f'{p}_rho'].ravel().astype(float)
    a2 = 4 * rho ** 2 * R[f'{p}_k'].ravel() * R[f'{p}_L'].ravel() ** 2 * S ** 2
    ok = np.isfinite(me) & np.isfinite(ae) & np.isfinite(S) & np.isfinite(G) & np.isfinite(a2)
    for q in cc.FIELDS:
        ok &= np.isfinite(c[q])
    z = lambda v: np.where(ok, v, 0.0)
    r = c['rho']
    d = dict(ok=ok, x=A, z=B, eta=mc.eta(R[f'{p}_T'].ravel()), cos=R[f'{p}_cos'].ravel(),
             nm=z(me * S ** 2), dm=z(S ** 2), na=z(ae * G ** 2), da=z(G ** 2), a2=z(a2),
             cnm=z(-r * c['uw'] * c['dUdz']), cdm=z(c['dUdz'] ** 2), cna=z(-r * c['wT'] * c['dTdz']), cda=z(c['dTdz'] ** 2),
             lnm=z(-r * c['uv'] * c['dUdy']), ldm=z(c['dUdy'] ** 2), lna=z(-r * c['vT'] * c['dTdy']), lda=z(c['dTdy'] ** 2),
             F2=z(r ** 2 * (c['uT'] ** 2 + c['vT'] ** 2 + c['wT'] ** 2)),
             g2=z(c['dTdx'] ** 2 + c['dTdy'] ** 2 + c['dTdz'] ** 2),
             Fg=z(-r * (c['uT'] * c['dTdx'] + c['vT'] * c['dTdy'] + c['wT'] * c['dTdz'])))
    with np.errstate(divide='ignore', invalid='ignore'):
        gg = c['dTdx'] ** 2 + c['dTdy'] ** 2 + c['dTdz'] ** 2
        d['cphi'] = -(c['uT'] * c['dTdx'] + c['vT'] * c['dTdy'] + c['wT'] * c['dTdz']) / np.sqrt(
            (c['uT'] ** 2 + c['vT'] ** 2 + c['wT'] ** 2) * gg)
        d['RwT'] = c['wT'] / (c['wrms'] * c['Trms'])
        d['RuT'] = c['uT'] / (c['urms'] * c['Trms'])
    d['bin'] = (np.floor(A + 2.0) if p == 'y0' else np.floor((A + 3.0) / 0.5)).astype(int).clip(0)
    return d


def sums(d, m):
    m = m & d['ok']
    nb = d['bin'].max() + 1
    return np.array([np.bincount(d['bin'][m], weights=d[k][m], minlength=nb) for k in KEYS]), int(m.sum())


def ratio(S, i):
    """Pr of the form (sum num_m / sum den_m) / (sum num_a / sum den_a) with rows i..i+3."""
    with np.errstate(divide='ignore', invalid='ignore'):
        return (S[i].sum() / S[i + 1].sum()) / (S[i + 2].sum() / S[i + 3].sum())


def boot(S, i, rng):
    have = np.nonzero(S[1] > 0)[0]
    if len(have) < 3:
        return np.nan, np.nan
    Sh = S[:, have]
    w = rng.multinomial(len(have), np.ones(len(have)) / len(have), size=NBOOT).astype(float)
    with np.errstate(divide='ignore', invalid='ignore'):
        pr = ((w @ Sh[i]) / (w @ Sh[i + 1])) / ((w @ Sh[i + 2]) / (w @ Sh[i + 3]))
    pr = pr[np.isfinite(pr)]
    if len(pr) < NBOOT // 2:
        return np.nan, np.nan
    return tuple(np.percentile(pr, [2.5, 97.5]))


def quart(v):
    v = v[np.isfinite(v)]
    return tuple(np.percentile(v, [25, 50, 75])) if v.size else (np.nan,) * 3


def cls(D6, D10, m6, m10, rng, lateral=False):
    S6, _ = sums(D6, m6)
    S10, n = sums(D10, m10)
    if n < 20:
        return None
    r = {'n': n}
    for b, S, d, m in ((6, S6, D6, m6), (10, S10, D10, m10)):
        mu = S[0].sum() / S[1].sum()
        r[f'pc{b}'] = ratio(S, 5)
        r[f'pt{b}'] = ratio(S, 0)
        r[f'pl{b}'] = ratio(S, 9)
        r[f'r2s{b}'] = 4 * mu ** 2 * S[1].sum() / S[4].sum() if S[4].sum() > 0 else np.nan
        N = S[15].sum()
        r[f'r2f{b}'] = N ** 2 / (S[14].sum() * S[13].sum()) if S[13].sum() > 0 else np.nan
        mm = m & d['ok']
        r[f'cphi{b}'] = quart(d['cphi'][mm])
        r[f'cth{b}'] = quart(d['cos'][mm])
    r['ci'] = boot(S10, 5, rng)
    r['cil'] = boot(S10, 9, rng)
    return r


def row(name, r, lat=False):
    if r is None:
        return f'{name:32s} (no points)'
    s = (f'{name:32s} n={r["n"]:6d} classic b006 {r["pc6"]:7.3f} b010 {r["pc10"]:7.3f} 95% [{r["ci"][0]:7.3f},{r["ci"][1]:7.3f}] '
         f'tensor b006 {r["pt6"]:7.3f} b010 {r["pt10"]:7.3f} | R2_flux {r["r2f6"]:.2f}/{r["r2f10"]:.2f} R2_stress {r["r2s6"]:.2f}/{r["r2s10"]:.2f} | '
         f'cos_phi b010 {r["cphi10"][1]:.2f}[{r["cphi10"][0]:.2f},{r["cphi10"][2]:.2f}] cos_theta b010 {r["cth10"][1]:.2f}[{r["cth10"][0]:.2f},{r["cth10"][2]:.2f}]')
    if lat:
        s += f' | lateral b006 {r["pl6"]:7.3f} b010 {r["pl10"]:7.3f} 95% [{r["cil"][0]:7.3f},{r["cil"][1]:7.3f}]'
    return s


def freestream_tu():
    """Tu = sqrt(2k/3)/U_inf, k = (u''^2 + v''^2 + w''^2)/2, at x/D -5 (tu read) and 0 (whole-plane read), z/D 2-3."""
    out = {}
    for lab, tag, x0, x1, shape in (('x/D=-5', 'tu', -5.1, -4.9, None), ('x/D=0', 'full', -0.1, 0.1, None)):
        f = os.path.join(cc.CORR_DIR, f'corr_{tag}_b010.npz')
        if not os.path.exists(f):
            continue
        z = np.load(f)
        P = dict(np.load(f'{mc.CACHE}/ref193_b010.npz')) if tag == 'full' else None
        if tag == 'full':
            x, zz = P['y0_A'].ravel()[z['idx']] / D, P['y0_B'].ravel()[z['idx']] / D
        else:
            x, zz = None, None
        k = z['k'].astype(float)
        if tag == 'tu':
            v = k[np.isfinite(k)]
        else:
            m = (x >= x0) & (x <= x1) & (zz >= 2) & (zz <= 3) & np.isfinite(k)
            v = k[m]
        tu = np.sqrt(2 * v / 3) / mc.U_INF
        log(f'Tu {lab}, z/D 2-3: median {np.median(tu) * 100:.2f} %, range {tu.min() * 100:.2f}-{tu.max() * 100:.2f} % (n {v.size})')


def main():
    rng = np.random.default_rng(1)
    R = {b: dict(np.load(f'{mc.CACHE}/ref193_b{b:03d}.npz')) for b in (6, 10)}
    Y = {b: prep(R[b], load_corr('full', b, R[b]['y0_A'].shape), 'y0') for b in R}
    X = {b: prep(R[b], load_corr('x10', b, R[b]['x10_A'].shape), 'x10') for b in R}
    res = {}

    def run(fam, name, f, which='y0', lat=False):
        A = Y if which == 'y0' else X
        res[(fam, name)] = cls(A[6], A[10], f(A[6]), f(A[10]), rng)
        log(row(f'{fam} {name}', res[(fam, name)], lat))

    xin = lambda x0, x1: (lambda d: (d['x'] >= x0) & (d['x'] < x1))
    freestream_tu()
    log('\n== Regions, whole z/D 0-3 and z/D >= 0.3 (centre plane); cross-plane x/D=10 (wall-normal and lateral) ==')
    for rn, x0, x1 in REGIONS:
        run('R', f'{rn} all z', lambda d, x0=x0, x1=x1: xin(x0, x1)(d))
        run('R', f'{rn} z>=0.3', lambda d, x0=x0, x1=x1: xin(x0, x1)(d) & (d['z'] >= 0.3))
    run('Rx', 'x10 all', lambda d: d['z'] >= 0, 'x10', True)
    run('Rx', 'x10 z>=0.3', lambda d: d['z'] >= 0.3, 'x10', True)
    log('\n== Z: wall-distance bands, centre plane ==')
    for rn, x0, x1 in REGIONS:
        for i, zl in enumerate(ZLAB):
            run('Z', f'{rn} z {zl}', lambda d, x0=x0, x1=x1, i=i: xin(x0, x1)(d) & (d['z'] >= ZB[i]) & (d['z'] < ZB[i + 1]))
    log('\n== Zx: wall-distance bands, cross-plane x/D=10 ==')
    for i, zl in enumerate(ZLAB):
        run('Zx', f'x10 z {zl}', lambda d, i=i: (d['z'] >= ZB[i]) & (d['z'] < ZB[i + 1]), 'x10', True)
    log('\n== E: eta bands (z/D >= 0.3) and wall layer, centre plane ==')
    for rn, x0, x1 in REGIONS:
        run('E', f'{rn} wall layer z<0.3', lambda d, x0=x0, x1=x1: xin(x0, x1)(d) & (d['z'] < 0.3))
        for en, lo, hi in EB:
            run('E', f'{rn} {en}', lambda d, x0=x0, x1=x1, lo=lo, hi=hi: xin(x0, x1)(d) & (d['z'] >= 0.3) & (d['eta'] >= lo) & (d['eta'] < hi))
    log('\n== Ex: eta bands, cross-plane x/D=10 ==')
    run('Ex', 'x10 wall layer z<0.3', lambda d: d['z'] < 0.3, 'x10', True)
    for en, lo, hi in EB:
        run('Ex', f'x10 {en}', lambda d, lo=lo, hi=hi: (d['z'] >= 0.3) & (d['eta'] >= lo) & (d['eta'] < hi), 'x10', True)
    log('\n== Classic Pr_t,R against tensor Pr_C in the film (methods appendix): z/D >= 0.3 ==')
    for rn, x0, x1 in REGIONS:
        r = res[('R', f'{rn} z>=0.3')]
        if r:
            log(f'{rn:10s} classic b010 {r["pc10"]:.3f}, tensor b010 {r["pt10"]:.3f}, ratio classic/tensor {r["pc10"] / r["pt10"]:.2f}')
    os.makedirs(cc.CORR_DIR, exist_ok=True)
    with open(TABLE, 'w') as fh:
        fh.write('\n'.join(LINES) + '\n')
    figure(res)


def pts(ax, x, r, col, ls, mk, dx, key='pc', grey=True):
    ok = [i for i, q in enumerate(r) if q is not None and np.isfinite(q[f'{key}10'])]
    xs = np.array([x[i] for i in ok], float)
    p10 = np.array([r[i][f'{key}10'] for i in ok])
    p6 = np.array([r[i][f'{key}6'] for i in ok])
    if key == 'pc':
        lo = np.array([r[i]['ci'][0] for i in ok])
        hi = np.array([r[i]['ci'][1] for i in ok])
        err = np.array([np.maximum(p10 - lo, 0), np.maximum(hi - p10, 0)])
        ax.errorbar(xs + dx, p10, yerr=np.where(np.isfinite(err), err, 0), color=col, ls=ls, marker=mk, capsize=2, zorder=3)
        ax.plot(xs + dx, p6, color=col, ls='none', marker=mk, mfc='none', zorder=3)
        if grey:
            ax.plot(xs + dx, [r[i]['pt10'] for i in ok], color='0.55', ls='none', marker=mk, ms=3, zorder=2)
    else:
        ax.plot(xs + dx, p10, color=col, ls=ls, marker=mk, zorder=3)
        ax.plot(xs + dx, p6, color=col, ls='none', marker=mk, mfc='none', zorder=3)


def figure(res):
    set_mpl_params(0.5, grid_on=False)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    col = {'jet': '#D55E00', 'near film': '#0072B2', 'far film': '#009E73', 'x/D = 10': '#000000'}
    sty = {'jet': (':', 'o', -0.15), 'near film': ('--', 's', -0.05), 'far film': ('-.', '^', 0.05), 'x/D = 10': ('-', 'D', 0.15)}
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.1))
    for k, ax in enumerate(axs):
        if k < 2:
            ax.axhspan(0.9, 1.2, color='0.85', zorder=0)
            ax.axhline(0.9, color='#000000', ls=':', zorder=1)
            ax.axhline(2.0, color='0.5', ls='-.', zorder=1)
            ax.set_ylim(-0.5, 3.0)
            ax.set_ylabel(r'$\mathrm{Pr}_{t,R}$')
        else:
            ax.set_ylim(-0.05, 1.05)
            ax.set_ylabel(r'$R^2$')
        ax.axhline(0.0, color='0.7', lw=0.5, zorder=0)
    ax = axs[0]
    for rn, _, _ in REGIONS:
        c, (ls, mk, dx) = col[rn], sty[rn]
        pts(ax, range(5), [res.get(('Z', f'{rn} z {zl}')) for zl in ZLAB], c, ls, mk, dx)
    pts(ax, range(5), [res.get(('Zx', f'x10 z {zl}')) for zl in ZLAB], col['x/D = 10'], '-', 'D', 0.15)
    ax.set_xticks(range(5))
    ax.set_xticklabels(ZLAB, rotation=40, ha='right')
    ax.set_xlabel(r'$z/D$ band')
    ax.text(0.96, 0.93, '(a)', transform=ax.transAxes, ha='right')
    ax = axs[1]
    names = ['wall layer z<0.3', 'core', 'mixing layer', 'edge']
    for rn, _, _ in REGIONS:
        c, (ls, mk, dx) = col[rn], sty[rn]
        pts(ax, range(4), [res.get(('E', f'{rn} {n}')) for n in names], c, ls, mk, dx)
    pts(ax, range(4), [res.get(('Ex', f'x10 {n}')) for n in names], col['x/D = 10'], '-', 'D', 0.15)
    ax.set_xticks(range(4))
    ax.set_xticklabels([r'wall, $z/D<0.3$', r'core, $\eta>0.4$', 'mixing, 0.1-0.4', 'edge, 0.02-0.1'], rotation=40, ha='right')
    ax.set_xlabel(r'$\eta$ band ($z/D\geq0.3$)')
    ax.text(0.96, 0.93, '(b)', transform=ax.transAxes, ha='right')
    ax = axs[2]
    for rn, _, _ in REGIONS:
        c, (ls, mk, dx) = col[rn], sty[rn]
        r = [res.get(('Z', f'{rn} z {zl}')) for zl in ZLAB]
        pts(ax, range(5), r, c, ls, mk, dx, key='r2f')
        pts(ax, range(5), r, c, '--', 'x', dx, key='r2s')
    ax.set_xticks(range(5))
    ax.set_xticklabels(ZLAB, rotation=40, ha='right')
    ax.set_xlabel(r'$z/D$ band')
    ax.text(0.96, 0.93, '(c)', transform=ax.transAxes, ha='right')
    h = [Line2D([], [], color=col['jet'], ls=':', marker='o', label='jet'),
         Line2D([], [], color=col['near film'], ls='--', marker='s', label='near film'),
         Line2D([], [], color=col['far film'], ls='-.', marker='^', label='far film'),
         Line2D([], [], color=col['x/D = 10'], ls='-', marker='D', label=r'cross-plane $x/D=10$'),
         Line2D([], [], color='k', ls='none', marker='o', label='b010 (10 FT)'),
         Line2D([], [], color='k', ls='none', marker='o', mfc='none', label='b006 (6 FT)'),
         Line2D([], [], color='0.55', ls='none', marker='o', ms=3, label='tensor form, b010'),
         Line2D([], [], color='k', ls=':', label='RANS 0.9'),
         Line2D([], [], color='0.5', ls='-.', label='K&B low FST, 2'),
         Line2D([], [], color='0.85', lw=6, label='K&B high FST, 0.9-1.2'),
         Line2D([], [], color='k', ls='-', marker='o', label=r'(c) $R^2_{\mathrm{flux}}$'),
         Line2D([], [], color='k', ls='--', marker='x', label=r'(c) $R^2_{\mathrm{stress}}$')]
    fig.legend(handles=h, loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.22, 1, 1))
    fig.savefig(OUT)
    plt.close(fig)
    log(f'saved {OUT}')


if __name__ == '__main__':
    main()
