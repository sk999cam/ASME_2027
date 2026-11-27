"""A posteriori tests: corrected RANS 50M runs against the 193M b010 (paper Fig. F20 and Table T2, Results 4.3).

Usage: python3 fig_aposteriori.py     (writes figures/fig_aposteriori.pdf and figures/code/data/aposteriori_table.txt)

Panels
  (a) laterally averaged eta-bar, (b) mid-pitch eta, both against x/D 0-20; (c) lateral profile of eta
  against |y/D|, plate points averaged over x/D 15-20. The centreline panel is not drawn (the table still lists it).
  Colour: 193M black, RANS 50M base blue (#0072B2); the other runs keep their colours.
  The curves use the functions and caches of fig_eta_x.py, so the 193M b010 and the baseline RANS curves are identical
  to that figure (centreline |y/D| < 0.1, mid-pitch ||y/D| - 2.91| < 0.1, bins of 0.1 D in x; eta-bar over one pitch).
  Panel (d): the provided STAR profile files (compare/*_eta_profile_15_20.npz) where the run is in one; otherwise the same
  definition from wall_plate.csv (plain mean of the plate faces per 0.1 D bin of |y/D|, YB of prt_sweep_eta.py). The 193M
  profile is computed from the plate-vertex cache of fig_eta_x (batch b010, same bins) and checked against the STAR file.

Adding a run: one entry in RUNS (key, folder, colour, line style, label, what it changes, profile file or None, in figure).
Runs not finished (A, gamma = 1, P1) are absent.
"""
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from functions import set_mpl_params  # noqa: E402
import fig_eta_x as FX  # noqa: E402

REPO = FX.REPO
STAR = REPO + '/starccm/rans_50M_Lm25b_BR15_FST_ZPG_p0dns'
OKABE_ITO = FX.OKABE_ITO
OUT = '/home/sk999/Documents/ASME_2027/figures/fig_aposteriori.pdf'
TABLE = '/home/sk999/Documents/ASME_2027/figures/code/data/aposteriori_table.txt'
XMAX = FX.XMAX
YB = np.linspace(0.0, 2.91, 30)
YC = 0.5 * (YB[1:] + YB[:-1])
QS = (5, 10, 15, 20)

# (key, folder suffix, colour, line style, label, what the run changes, profile npz suffix/label or None, in figure)
RUNS = [
    ('base', '', '#0072B2', '-', r'RANS 50M, SST, $\Pr_t$ 0.9 (base)', 'none (baseline)',
     ('substB', 'SST, Pr_t 0.9'), True),
    ('pr04', '_prt04', '#56B4E9', '--', r'$\Pr_t$ 0.4', 'Pr_t 0.4, whole domain',
     ('substB', 'SST, Pr_t 0.4'), True),
    ('pr05', '_prt05', '#F0E442', '--', r'$\Pr_t$ 0.5', 'Pr_t 0.5, whole domain',
     ('substB', 'SST, Pr_t 0.5'), False),
    ('pr03', '_prt03', '#F0E442', ':', r'$\Pr_t$ 0.3', 'Pr_t 0.3, whole domain',
     ('substB', 'SST, Pr_t 0.3'), False),
    ('qcr', '_qcr', '#D55E00', '-.', 'QCR', 'stress relation (QCR 2000)',
     ('substB', 'SST + QCR, Pr_t 0.9'), False),
    ('qcr04', '_qcr_prt04', '#CC79A7', ':', r'QCR, $\Pr_t$ 0.4', 'QCR and Pr_t 0.4, whole domain', None, True),
    ('C', '_substC', '#D55E00', '--', 'run C (3DNS heat flux)',
     '50M 3DNS heat flux imposed, all regions',
     ('substB', 'SST + 3DNS heat flux (run C)'), False),
    ('B', '_substB', '#0072B2', '-', 'run B (3DNS stresses and heat flux)',
     '50M 3DNS stresses and heat flux imposed, all regions',
     ('substB', 'SST + 3DNS stresses and heat flux (run B)'), False),
    ('P2', '_substCfilm', '#009E73', (0, (3, 1, 1, 1, 1, 1)), 'resolved heat flux in the film', '50M 3DNS heat flux imposed, x/D 3-25 only',   # colour/style not #E69F00 '-.': that is the 20M LES (paper-wide)
     ('substCfilm', 'SST + 3DNS heat flux x/D 3-25 (P2)'), True),
    # prepared, inactive until the run folder holds post_run1/wall_plate.csv and convergence/convergence.md
    ('gamma1', '_gamma1', '#D55E00', '-', r'$\gamma=1$', 'intermittency gamma = 1', None, True),
    ('P1', '_prtfilm', '#0072B2', (0, (3, 1, 1, 1, 1, 1)), r'P1 ($\Pr_t$ 0.4, film)', 'Pr_t 0.4 in the film only', None, True),
]
NOTES = {'B': 'imposes 50M 3DNS inputs (C50b b006); BR 1.416; T up to 303.6 K outside x/D 0-20; approach layer changed',
         'C': 'imposes 50M 3DNS inputs (C50b b006)',
         'P2': 'imposes 50M 3DNS inputs (C50b b006)'}
RUNS = [r for r in RUNS if os.path.exists(STAR + r[1] + '/post_run1/wall_plate.csv')
        and os.path.exists(STAR + r[1] + '/convergence/convergence.md')]
REF_NAME = '193M b010'


def folder(r):
    return STAR + r[1]


def plate(r):
    w = FX.read_csv(folder(r) + '/post_run1/wall_plate.csv')
    return w['X'] / FX.D, w['Y'] / FX.D, FX.eta_of(w['Temperature'])


def profile_plain(x, y, e):
    s = (x >= 15) & (x <= 20) & np.isfinite(e)
    k = np.digitize(np.abs(y[s]), YB) - 1
    return np.array([np.nanmean(e[s][k == i]) if np.any(k == i) else np.nan for i in range(len(YB) - 1)])


def profile_file(spec):
    suf, lab = spec
    z = np.load(f'{STAR}_{suf}/compare/{suf}_eta_profile_15_20.npz', allow_pickle=True)
    i = list(z['labels']).index(lab)
    return z['y_over_D'], z['eta'][i]


def convergence(r):
    t = open(folder(r) + '/convergence/convergence.md').read()
    it = int(re.search(r'at iteration (\d+)', t).group(1))
    res = [(m.group(1), float(m.group(2))) for m in re.finditer(r'^\| (\S+) \| ([0-9.e+-]+) \| [0-9.]+ \| [+-][0-9.]+ \|', t, re.M)]
    br = float(re.search(r'^\| BR \| ([0-9.]+)', t, re.M).group(1))
    mx = max(res, key=lambda p: p[1])
    return it, mx, br


def quantities(xb, lat, xc, cl, mp):
    d = {'bar': FX.window_means(xb, lat)}
    d['cl'] = [FX.at(xc, cl, q) for q in QS]
    d['mp'] = [FX.at(xc, mp, q) for q in QS]
    return d


def main():
    # Drawn at the printed size for full text width (7.5 in; figure* with width=\textwidth, scale 1); body text 10 pt.
    set_mpl_params(10 / 22, grid_on=True)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    fig, axs = plt.subplots(1, 3, figsize=(7.5, 2.4))   # centreline panel dropped (audit 2026-11-27); the table still holds the centreline values
    ax_bar, ax_mp, ax_pr = axs
    data = {}

    # reference 193M b010 (identical to fig_eta_x)
    x, y, e = FX.plate_points('193M')
    xc, cl, mp = FX.lines(x, y, e)
    xb, lat = FX.lateral('193M')
    ref = quantities(xb, lat, xc, cl, mp)
    ref_prof = profile_plain(x, y, e)
    zf = np.load(STAR + '_substB/compare/substB_eta_profile_15_20.npz', allow_pickle=True)
    ref_file = zf['eta'][list(zf['labels']).index('193M 3DNS (C193b b010)')]
    print('193M profile, max |computed - STAR file|: %.2e' % np.nanmax(np.abs(ref_prof - ref_file)))

    def draw(xb_, lat_, xc_, cl_, mp_, prof, **kw):
        s = np.isfinite(lat_) & (xb_ > 2.36) & (xb_ <= XMAX)
        ax_bar.plot(xb_[s], lat_[s], **kw)
        m = np.isfinite(mp_) & (xc_ > 2.36) & (xc_ <= XMAX)
        ax_mp.plot(xc_[m], mp_[m], **kw)
        ax_pr.plot(YC, prof, **kw)

    # RANS runs
    for r in RUNS:
        key, _, col, ls, lab, _, pspec, infig = r
        x, y, e = plate(r)
        xc, cl, mp = FX.lines(x, y, e)
        xb, lat = FX.lateral_mean(x * FX.D, y * FX.D, e, np.arange(-2 * FX.D, 25.5 * FX.D, 0.1 * FX.D))
        xb = xb / FX.D
        prof_plain = profile_plain(x, y, e)
        if pspec:
            yy, prof = profile_file(pspec)
            note = 'STAR file; max |file - plain mean| %.1e' % np.nanmax(np.abs(prof - prof_plain))
        else:
            prof, note = prof_plain, 'computed from wall_plate.csv (no STAR profile file)'
        print(f'{key:6s} profile: {note}')
        data[key] = dict(quantities(xb, lat, xc, cl, mp), prof=prof, conv=convergence(r), note=note)
        if infig:
            draw(xb, lat, xc, cl, mp, prof, color=col, ls=ls, label=lab)

    # 20M LES (resolution ladder, owner 2026-11-27): same definitions and caches as fig_eta_x ('20M LES', batch b006, plate-vertex cache);
    # style of fig_eta_x ('#E69F00', dash-dot), which is also the style of run P2 in this figure (clash, kept as ruled).
    LES_KEY = '20M LES'
    x, y, e = FX.plate_points(LES_KEY)
    xc, cl, mp = FX.lines(x, y, e)
    xb, lat = FX.lateral(LES_KEY)
    les_prof = profile_plain(x, y, e)
    data['LES'] = dict(quantities(xb, lat, xc, cl, mp), prof=les_prof)
    lc, ls_ = FX.STYLE[LES_KEY]
    draw(xb, lat, xc, cl, mp, les_prof, color=lc, ls=ls_, label=LES_KEY)

    # reference drawn last, thick black
    xb, lat = FX.lateral('193M')
    xc, cl, mp = FX.lines(*FX.plate_points('193M'))
    kw = dict(color='k', ls='-', lw=2.5 * plt.rcParams['lines.linewidth'], label='193M', zorder=1)
    draw(xb, lat, xc, cl, mp, ref_prof, **kw)

    # experiment id727, as fig_eta_x draws it (markers only); lateral profile from the same map, z folded onto |y|
    xe = np.loadtxt(FX.EXP + 'x_D.txt')[0]
    ee = np.nanmean(np.loadtxt(FX.EXP + FX.EXP_FILE), axis=0)
    ex, ey, ev = FX.exp_points()
    exl = FX.lines(ex, ey, ev)
    mk = dict(marker='o', ls='none', color='k', mfc='none', ms=3.5, label='experiment')
    s = (xe >= 2.2) & (xe <= XMAX) & np.isfinite(ee)
    ax_bar.plot(xe[s][::8], ee[s][::8], **mk)
    for ax, i in ((ax_mp, 2),):
        m = np.isfinite(exl[i]) & (exl[0] >= 2.2) & (exl[0] <= XMAX)
        ax.plot(exl[0][m][::4], exl[i][m][::4], **mk)
    exp_prof = profile_plain(ex, ey, ev)
    ax_pr.plot(YC, exp_prof, **mk)
    exp = {'bar': FX.window_means(xe, ee), 'cl': [FX.at(exl[0], exl[1], q) for q in QS],
           'mp': [FX.at(exl[0], exl[2], q) for q in QS]}

    ax_bar.set_ylim(0, 0.45)
    ax_pr.set_ylim(bottom=0)
    ax_mp.set_ylim(0, 0.06)
    for ax in (ax_bar, ax_mp):
        ax.set_xlim(0, XMAX)
        ax.set_xlabel(r'$x/D$')
    ax_pr.set_xlim(0, 2.91)
    ax_pr.set_xlabel(r'$|y/D|$')
    for ax, lab, yl in zip(axs, '(a) (b) (c)'.split(),
                           (r'$\bar{\eta}$', r'$\eta$ (mid-pitch)', r'$\eta$ ($x/D$ 15--20)')):
        ax.set_ylabel(yl)
        ax.text(0.0, 1.02, lab, transform=ax.transAxes, va='bottom', ha='left')
    h, l = ax_bar.get_legend_handles_labels()
    o = sorted(range(len(l)), key=lambda i: (l[i] != '193M', l[i] != 'experiment', i))
    h, l = [h[i] for i in o], [l[i] for i in o]
    fig.tight_layout(pad=0.3, rect=(0, 0.2, 1, 1))
    fig.legend(h, l, loc='lower center', ncol=3, frameon=False, handlelength=2.0, columnspacing=0.8, labelspacing=0.25,
               fontsize=plt.rcParams['font.size'] - 1, bbox_to_anchor=(0.5, 0.0))
    fig.savefig(OUT)
    fig.savefig('/home/sk999/Documents/ASME_2027/planning/review/55_fig_aposteriori_les.png', dpi=200)

    write_table(ref, ref_prof, data, exp)


def write_table(ref, ref_prof, data, exp):
    f = lambda v: f'{v:.3f}'
    L = ['Table T2: a posteriori runs, RANS 50M against the 193M b010 (fig_aposteriori.py)',
         'eta = (T0,fs - T_w)/(T0,fs - T0,c). eta-bar: spanwise-weighted lateral mean (bl_post.lateral_mean), mean of the 0.1 D bins',
         'over x/D 2.5-7, 7-15, 15-20. Centreline |y/D| < 0.1, mid-pitch ||y/D| - 2.91| < 0.1, bins of 0.1 D in x; value at x/D 5, 10, 15',
         '= mean of the bins within +-0.26 D, at x/D 20 = bins 19.5-20.0 (fig_eta_x.at). Nothing beyond x/D 20.',
         'Difference = run - 193M b010. Iteration and residual from convergence/convergence.md (largest normalised residual of the', 'nine equations at that iteration); BR from the same file.', '']
    names = {r[0]: r for r in RUNS}
    hdr = ('run', 'bar2.5-7', 'bar7-15', 'bar15-20', 'cl5', 'cl10', 'cl15', 'cl20', 'mp5', 'mp10', 'mp15', 'mp20')
    L.append('Values (eta)')
    L.append(' '.join(f'{h:>9s}' for h in hdr))
    def row(name, d):
        return f'{name:9s} ' + ' '.join(f'{v:9.3f}' for v in d['bar'] + d['cl'] + d['mp'])
    L.append(row(REF_NAME[:9], ref))
    L.append(row('exp id727', exp) + '   (7-7-7, DR 1.5, BR 1.5, Tu 0.5 %; fig_eta_x definitions; no iteration/residual)')
    for r in RUNS:
        L.append(row(r[0], data[r[0]]))
    L.append(row('LES 20M', data['LES']) + '   (20M LES C20bW stat3 b006, fig_eta_x definitions and plate-vertex cache; no iteration/residual)')
    L += ['', 'Differences from the 193M b010 (run - 193M)']
    L.append(' '.join(f'{h:>9s}' for h in hdr))
    for nm, d in [(r[0], data[r[0]]) for r in RUNS] + [('LES 20M', data['LES'])]:
        v = [a - b for a, b in zip(d['bar'] + d['cl'] + d['mp'], ref['bar'] + ref['cl'] + ref['mp'])]
        L.append(f'{nm:9s} ' + ' '.join(f'{x:+9.3f}' for x in v))
    pr = lambda p: np.nanmax(np.abs(p - ref_prof))     # noqa: E731
    L += ['', 'Lateral profile x/D 15-20 (panel c), max |run - 193M| over the 29 bins of |y/D| (eta):']
    L += [f'  {r[0]:9s} {pr(data[r[0]]["prof"]):.3f}' for r in RUNS] + [f'  {"LES 20M":9s} {pr(data["LES"]["prof"]):.3f}']
    L += ['', 'Convergence and what each run changes',
          f'{"run":9s} {"iteration":>9s} {"max resid.":>10s} {"equation":>13s} {"BR":>7s}  what the run changes']
    for r in RUNS:
        it, (eq, rv), br = data[r[0]]['conv']
        L.append(f'{r[0]:9s} {it:9d} {rv:10.2e} {eq:>13s} {br:7.4f}  {r[5]}  | notes: {NOTES.get(r[0], "-")}')
    L += ['', 'Run keys: base = rans_50M_Lm25b_BR15_FST_ZPG_p0dns; pr04/pr05/pr03 = _prt04/_prt05/_prt03; qcr = _qcr; qcr04 = _qcr_prt04;',
          'B = _substB (3DNS stresses and heat flux, 50M C50b b006); P2 = _substCfilm (3DNS heat flux in regions 5 and 6, x/D 3-25).',
          'C = _substC (3DNS heat flux, all regions). Runs A, gamma = 1 and P1 are not finished and are absent; the figure draws base, pr04, qcr04, P2 only.', '',
          'Flags']
    L += FLAGS
    os.makedirs(os.path.dirname(TABLE), exist_ok=True)
    with open(TABLE, 'w') as fh:
        fh.write('\n'.join(L) + '\n')
    print('\n'.join(L))


FLAGS = [
    '1. Runs B and P2 impose stresses / heat flux from the 50M 3DNS (C50b b006), not the 193M; all differences are against the 193M b010.',
    '2. STAR tables (substB_eta_windows_x20.txt etc.) use mid-pitch ||y/D| - 2.91| < 0.15 and centreline |y/D| < 0.1, x/D 5 and 10 windows',
    '   +-0.25; this table follows fig_eta_x (mid-pitch +-0.1, x/D 5, 10, 15 windows +-0.26). Mid-pitch values differ by a few 1e-4 to 1e-3.',
    '3. qcr_prt04/compare/qcr_prt04_eta_windows.txt is an older table: its x/D 20 window differs (baseline centreline 0.4936 against',
    '   0.4969 in the x/D 19.5-20.0 tables). It is not used here; the QCR + Pr_t 0.4 values are recomputed from wall_plate.csv.',
    '4. QCR + Pr_t 0.4 has no STAR profile file; its profile is recomputed from wall_plate.csv with the STAR definition.',
    '5. Run B has BR 1.416 against 1.527 (193M) and 1.59 for the other runs (the stresses cut the discharge at fixed p0).',
    '6. Run B has T up to 303.6 K (cells above T0 at the wall, none in x/D 0-20); its approach layer differs (Re_theta 224 at x/D -2).',
]

if __name__ == '__main__':
    main()
