#!/usr/bin/env python3
"""F9 fig_spots: turbulent spots between the jets, 193M 3DNS (b010, left column) and 20M LES (b006, right column).

Detector (one definition for every intermittency in the paper, helper _gamma_extract_v.py): d = (v - <v>_T)^2 / U_inf^2 at
z/D = 0.1, v the spanwise velocity, <v>_T the running mean in time over T = 2 D/U, d smoothed by a 1 D/U running mean;
laminar / turbulent by Otsu's method on log10 d (pooled over the intermittent band |y/D| 2 to 3, x/D 0 to 20, of the 193M; the
same threshold is applied to both runs); gamma = fraction of time above the threshold. Definition and threshold policy: planning/maths_audit/intermittency.
Top, four plan views of log10 d (z/D 0.1, x/D -3 to 20, y/D folded so that the laminar wedge between two jets is centred at
y/D 3): white contour = threshold, hole exit outlined (with its periodic copy), dashed line = mid-pitch line. The times follow one
turbulent patch (the longest-lived spot of the 40 D/U window) as it moves downstream.
Middle, (x, t) diagram of log10 d on the mid-pitch line y/D = 3 over 40 D/U; contour = threshold; ticks at the right edge mark the
snapshot times t1..t4. Bottom, gamma(x) on the mid-pitch line (time fraction above the threshold over the batch window), both runs.
One colour scale (inferno) for all maps. usage: python3 fig_spots.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402
from cycler import cycler  # noqa: E402

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from functions import set_mpl_params  # noqa: E402
from hole_exit_outline import load as load_hole  # noqa: E402
import fig_wall_state as WS  # noqa: E402      # RANS 50M loader (rans_slab, regrid)
import _gamma_extract as GX  # noqa: E402      # old detector, printed comparison only
import _gamma_extract_v as GV  # noqa: E402
from scipy.ndimage import label  # noqa: E402

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_spots.png'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
XL = (-3.0, 20.0)
XH = (-1.0, 2.35)
CASES = [('c193b', '193M 3DNS'), ('c20bw', '20M LES')]
TWIN = 40.0                    # length of the (x, t) window, D/U
NSNAP = 4
FW = 7.5                       # printed width: figure* at \textwidth, scale 1 (height computed in main); font 10 pt = body size
BG = dict(facecolor='k', edgecolor='none', pad=1.5, alpha=0.75)
XQ = (-2, 5, 10, 15, 20)


def find_spots(mask, t, x, minrec=30):
    """connected regions of the binary mid-pitch (x, t) field that are born narrow (<= 0.6 D, centre x/D 0 to 14), are followed
    while narrower than 5 D and last >= 6 D/U; dicts with t0, x0, life and the edge trajectories (t, le, te). Used to choose
    the patch that is marked in the plan views; no celerity is derived (see fig_spots_stats_v.py, audit only)."""
    lab, n = label(mask, structure=np.ones((3, 3), int))
    out = []
    for k in range(1, n + 1):
        rr, ii = np.nonzero(lab == k)
        rows = np.unique(rr)
        if len(rows) < minrec:
            continue
        le = np.array([x[ii[rr == r]].max() for r in rows])
        te = np.array([x[ii[rr == r]].min() for r in rows])
        if (le[0] - te[0]) > 0.6 or not (0.0 <= 0.5 * (le[0] + te[0]) <= 14.0):
            continue
        bad = np.where((le - te) >= 5.0)[0]
        stop = bad[0] if len(bad) else len(rows)
        if stop < minrec or t[rows[stop - 1]] - t[rows[0]] < 6.0:
            continue
        out.append(dict(t0=float(t[rows[0]]), x0=float(te[0]), life=float(t[rows[stop - 1]] - t[rows[0]]), t=t[rows[:stop]],
                        le=le[:stop], te=te[:stop]))
    return out


def fold_plane(y, f):
    """f (nx, ny) on y -3..3 -> yf 0..6 (y < 0 shifted by 6), sorted; the periodic duplicate y = 3 stays at 3."""
    yf = np.concatenate([y[y >= 0], y[y < 0] + 6.0])
    ff = np.concatenate([f[:, y >= 0], f[:, y < 0]], axis=1)
    o = np.argsort(yf, kind='stable')
    return yf[o], ff[:, o]


def pick_window(tag, lt):
    """the spot of the batch with the largest downstream travel (x of its leading edge at the end minus its birth x, x/D <= 20);
    returns (window start, spot dict, snapshot times, cache dict, line field)."""
    Z = GV.load(tag, 'z010')
    t, ok = Z['tu'], GV.XG <= 20.0 + 1e-9
    ld = Z['line_y+3.00'][:, ok]
    sp = find_spots(ld > lt, t, GV.XG[ok])
    p0, p1 = t[0] + 3.0, t[-1] - 3.0
    sp = [s for s in sp if s['t0'] >= p0 and s['t'][-1] <= p1]
    nmax = int((TWIN - 6.0) / GV.DT)                      # follow the spot for at most 34 D/U
    s = max(sp, key=lambda s: s['le'][min(len(s['le']) - 1, nmax)] - s['x0'])
    w0 = float(np.clip(s['t0'] - 4.0, p0, p1 - TWIN))
    life = min(s['life'], TWIN - 6.0)
    tk = s['t0'] + life * np.array([0.08, 0.38, 0.68, 0.97])[:NSNAP]
    return w0, s, tk, Z, ld


def main():
    set_mpl_params(10.0 / 22.0, grid_on=False)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    HX, HY = load_hole()
    D = {}
    shown = []
    for tag, name in CASES:
        lt = GV.policy_logth(tag)
        w0, sp, tk, Z, ld = pick_window(tag, lt)
        tf, frames = GV.plane_logd(tag, w0, w0 + TWIN)
        ik = [int(np.argmin(np.abs(tf - q))) for q in tk]
        tt = Z['tu']
        m = (tt >= w0) & (tt <= w0 + TWIN)
        D[tag] = dict(lt=lt, w0=w0, sp=sp, tk=tf[ik], ik=ik, frames=frames, tt=tt[m], line=ld[m], Z=Z)
        shown += [np.asarray(frames[i], np.float32).ravel()[::7] for i in ik] + [ld[m].ravel()[::3]]
        assert tf[ik].min() >= w0 and tf[ik].max() <= w0 + TWIN, 'snapshot times outside the (x, t) window'
        print(f'{tag}: threshold log10 d {lt:.3f}; window t U/D {w0:.1f}-{w0 + TWIN:.1f}; tracked patch born t {sp["t0"]:.1f} '
              f'x/D {sp["x0"]:.2f}, life {sp["life"]:.1f} D/U; snapshot times {np.round(tf[ik], 1)}', flush=True)
    vmin, vmax = np.percentile(np.concatenate(shown), [1, 99.5])
    print(f'colour range log10 d: {vmin:.2f} to {vmax:.2f}')

    FH = 0.08 + NSNAP * 0.62 + (NSNAP - 1) * 0.14 + 0.45 + 1.8 + 0.62 + 1.35 + 0.55    # inch
    fig = plt.figure(figsize=(FW, FH))
    LM, RM, GAP = 0.62 / FW, 1 - 0.85 / FW, 0.2 / FW
    W = (RM - LM - GAP) / 2
    cols = [LM, LM + W + GAP]
    inch = lambda a: a / FH
    PH, PG = inch(0.62), inch(0.14)
    H_XT, H_G = inch(1.8), inch(1.35)
    top0 = 1 - inch(0.08)
    p_bot = top0 - NSNAP * PH - (NSNAP - 1) * PG
    xt_top = p_bot - inch(0.45)
    xt_bot = xt_top - H_XT
    g_top = xt_bot - inch(0.62)
    g_bot = g_top - H_G
    im = None
    for ci, (tag, name) in enumerate(CASES):
        d = D[tag]
        lt = d['lt']
        x = GV.XG
        sp = d['sp']
        for k, ik in enumerate(d['ik']):
            top = top0 - k * (PH + PG)
            ap = fig.add_axes([cols[ci], top - PH, W, PH])
            yf, ff = fold_plane(GV.YG, np.asarray(d['frames'][ik], np.float32))
            im = ap.pcolormesh(x, yf, ff.T, cmap='inferno', vmin=vmin, vmax=vmax, shading='nearest', rasterized=True)
            ap.contour(x, yf, ff.T, levels=[lt], colors='w', linewidths=0.8, linestyles='solid', zorder=3)
            for dy in (0.0, 6.0):
                ap.add_patch(Polygon(np.c_[HX, HY + dy], closed=True, facecolor='w', edgecolor='k', lw=0.8, zorder=4))
            ap.plot([XL[0], XL[1]], [3.0, 3.0], color=OKABE_ITO[2], ls='--', lw=0.9, zorder=5)
            ap.set_xlim(*XL)
            ap.set_ylim(0, 6)
            ap.set_yticks([0, 3, 6])
            ap.text(0.012, 0.9, r'$t_%d=%.1f$' % (k + 1, d['tk'][k] - d['w0']), transform=ap.transAxes, color='w',
                    ha='left', va='top', zorder=6, bbox=BG)
            if k == 0:
                ap.text(0.99, 0.9, '(%s) %s' % ('ab'[ci], name), transform=ap.transAxes, color='w', ha='right', va='top',
                        zorder=6, bbox=BG)
            if ci == 0:
                ap.set_ylabel(r'$y/D$', labelpad=1)
            else:
                ap.set_yticklabels([])
            if k < NSNAP - 1:
                ap.set_xticklabels([])
            else:
                ap.set_xlabel(r'$x/D$', labelpad=1)
            tcur = d['tk'][k]
            j = int(np.argmin(np.abs(sp['t'] - tcur)))
            xc = 0.5 * (sp['le'][j] + sp['te'][j])             # marker of the tracked patch (centre of its x-extent)
            ap.plot([xc], [2.35], marker='^', ms=5, mfc='w', mec='k', mew=0.7, ls='none', zorder=8)
            if k == 0 and ci == 1:
                ap.annotate('tracked patch', xy=(xc, 2.25), xytext=(min(xc + 4.5, 17.0), 1.0), color='w',
                            ha='center', va='center', bbox=BG, arrowprops=dict(arrowstyle='->', color='w', lw=0.8), zorder=7)
            if k == 0 and ci == 0:
                ap.annotate('tracked patch', xy=(xc, 2.25), xytext=(min(xc + 4.5, 17.0), 1.0), color='w',
                            ha='center', va='center', bbox=BG, arrowprops=dict(arrowstyle='->', color='w', lw=0.8), zorder=7)
            if ci == 0 and k == 0:
                ap.annotate('mid-pitch', xy=(4.9, 3.0), xytext=(4.9, 1.0), color='w', ha='center',
                            va='center', bbox=BG, arrowprops=dict(arrowstyle='->', color='w', lw=0.8), zorder=7)
        # (x, t) diagram
        ax = fig.add_axes([cols[ci], xt_bot, W, H_XT])
        xl = x[:d['line'].shape[1]]                            # the line field stops at x/D 20
        tt = d['tt'] - d['w0']                                  # time relative to the window start t0 of each run
        dtm = tt[1] - tt[0]
        ax.imshow(d['line'], origin='lower', aspect='auto', cmap='inferno', vmin=vmin, vmax=vmax, interpolation='nearest',
                  extent=[xl[0] - 0.01, xl[-1] + 0.01, tt[0] - dtm / 2, tt[-1] + dtm / 2])
        ax.contour(xl, tt, d['line'], levels=[lt], colors='w', linewidths=0.8, linestyles='solid')
        for xh in XH:
            ax.axvline(xh, color='w', ls='--', lw=0.9)
        ax.set_xlim(*XL)
        ax.set_ylim(tt[0], tt[-1])
        ax.text(0.5 * sum(XH), tt[-1] - 0.03 * TWIN, 'hole\nexit', color='w', ha='center', va='top')
        for k, tk in enumerate(d['tk'] - d['w0']):             # snapshot times (relative): ticks on the left edge, one label column
            ax.plot([XL[0], XL[0] + 0.7], [tk, tk], color='w', lw=1.2, zorder=5)
            if k in (0, NSNAP - 1):                            # the ticks in between are too close for 10 pt labels
                ax.text(XL[0] + 0.9, tk, r'$t_%d$' % (k + 1), color='w', ha='left', va='top' if k == 0 else 'bottom',
                        bbox=dict(facecolor='k', edgecolor='none', pad=0.5, alpha=0.85), zorder=6)
        ax.text(0.99, 0.97, '(%s) %s' % ('cd'[ci], name), transform=ax.transAxes, color='w', ha='right', va='top',
                bbox=BG)
        ax.text(0.99, 0.03, r'$t_0U_\infty/D$ = %.1f' % d['w0'], transform=ax.transAxes, color='w', ha='right', va='bottom',
                bbox=BG)
        ax.set_xlabel(r'$x/D$', labelpad=1)
        if ci == 0:
            ax.set_ylabel(r'$(t-t_0)U_\infty/D$', labelpad=1)
        else:
            ax.set_yticklabels([])
    cax = fig.add_axes([RM + 0.012, xt_bot, 0.014, top0 - xt_bot])
    fig.colorbar(im, cax=cax).set_label(r'$\log_{10}d$, $d=(v-\langle v\rangle_T)^2/U_\infty^2$ ($z/D=0.1$)')
    # gamma(x) on the mid-pitch line, both runs
    ag = fig.add_axes([LM, g_bot, RM - LM, H_G])
    sty = {'c193b': (OKABE_ITO[5], '-'), 'c20bw': (OKABE_ITO[6], '--')}
    print('mid-pitch gamma (line y/D = 3) at x/D ' + ' '.join(map(str, XQ)))
    for tag, name in CASES:
        d = D[tag]
        g = GV.load_z010(tag)['gamma']
        gl = 0.5 * (g[:, 0] + g[:, -1])
        ok = GV.XG <= XL[1] + 1e-9
        fr = GV.load(tag, 'z010')['frac'][1]
        gb = [GV.gamma_at(fr, d['lt'] + np.log10(f)) for f in GV.FAC]       # threshold x0.5 (upper) and x2 (lower)
        gb = [0.5 * (a[:, 0] + a[:, -1]) for a in gb]
        ag.fill_between(GV.XG[ok], gb[1][ok], gb[0][ok], color=sty[tag][0], alpha=0.15, lw=0)
        for ge in gb:                                     # thin edge lines in the line colour keep the two bands apart
            ag.plot(GV.XG[ok], ge[ok], color=sty[tag][0], lw=0.5, alpha=0.8, ls=sty[tag][1])
        if tag == 'c193b':
            band_xy = (12.5, float(np.interp(12.5, GV.XG, gb[1])))     # lower (threshold x2) edge of the 193M band
        print(f'  {name}: threshold x0.5 / x2 at x/D {XQ}: ' + ' | '.join(
            f'{np.interp(q, GV.XG, gb[0]):.3f} / {np.interp(q, GV.XG, gb[1]):.3f}' for q in XQ))
        ag.plot(GV.XG[ok], gl[ok], color=sty[tag][0], ls=sty[tag][1], lw=1.6)
        go = GX.load_z010(tag)['gamma'][0][:, 300]
        gli = GV.gamma_at(fr, GV.policy_logth(tag, 2.0, 'iii'))
        gli = 0.5 * (gli[:, 0] + gli[:, -1])
        print(f'  {name}: old |du/dt| ' + ' '.join(f'{np.interp(q, GX.XG, go):.3f}' for q in XQ)
              + ' -> line policy ' + ' '.join(f'{np.interp(q, GV.XG, gli):.3f}' for q in XQ)
              + ' -> band policy ' + ' '.join(f'{np.interp(q, GV.XG, gl):.3f}' for q in XQ))
        for nm_, a_, o_ in (('old', go, GX.XG), ('line', gli, GV.XG), ('band', gl, GV.XG)):
            m_ = (o_ >= 0) & (o_ <= 20 + 1e-9) & (a_[:len(o_)] > 0.8)
            print(f'     first x/D with mid-pitch gamma > 0.8 ({nm_}): ' + (f'{o_[m_][0]:.2f}' if m_.any() else 'none'))
        if tag == 'c193b':
            i = np.argmin(np.abs(GV.XG - 14.0))
            ag.annotate('193M 3DNS', xy=(GV.XG[i], gl[i]), xytext=(12.2, 0.93), color=sty[tag][0], ha='center',
                        va='center', arrowprops=dict(arrowstyle='->', color=sty[tag][0], lw=0.8))
        else:
            i = np.argmin(np.abs(GV.XG - 10.5))
            ag.annotate('20M LES', xy=(GV.XG[i], gl[i]), xytext=(6.6, 0.4), color=sty[tag][0], ha='center',
                        va='center', arrowprops=dict(arrowstyle='->', color=sty[tag][0], lw=0.8))
    # RANS 50M: transported intermittency (gamma-Re_theta model variable), z/D 0.1 layer, mid-pitch line y/D = 3; the same
    # loader and cache as fig_wall_state.py panel (d). A model variable, not the detector gamma: no threshold band.
    r = WS.rans_slab()
    xr = np.arange(XL[0], XL[1] + 1e-9, 0.04)
    sel = (r['z'] > 0.094) & (r['z'] < 0.102)
    gr = WS.regrid(r['x'][sel], r['y'][sel], r['g'][sel], xr, np.array([3.0]), 0.15)[:, 0]
    ag.plot(xr, gr, color=OKABE_ITO[3], ls='-.', lw=1.6)
    i = np.nanargmin(np.abs(xr - 15.0) + np.where(np.isfinite(gr), 0, 1e3))
    ag.annotate('RANS 50M (model $\\gamma$)', xy=(xr[i], gr[i]), xytext=(16.0, 0.42), color=OKABE_ITO[3],
                ha='center', va='center', arrowprops=dict(arrowstyle='->', color=OKABE_ITO[3], lw=0.8))
    print('RANS 50M model gamma, mid-pitch y/D 3, z/D 0.1 layer: x/D 5 10 15 20: '
          + ' '.join(f'{np.interp(q, xr[np.isfinite(gr)], gr[np.isfinite(gr)]):.3f}' for q in (5, 10, 15, 20)))
    print('RANS model gamma at mid-pitch x/D 15 by cached layer (z/D 0.05 and 0.2 are not cached; cells only 0.08 < z/D < 0.12):')
    for lo_, hi_ in ((0.080, 0.087), (0.087, 0.094), (0.094, 0.102), (0.102, 0.109), (0.109, 0.12)):
        s_ = (r['z'] > lo_) & (r['z'] <= hi_)
        gg_ = WS.regrid(r['x'][s_], r['y'][s_], r['g'][s_], np.array([15.0]), np.array([3.0]), 0.15)[0, 0]
        print(f'   z/D {r["z"][s_].mean():.3f}: {gg_:.3f}')
    for xh in XH:
        ag.axvline(xh, color='0.4', ls=':', lw=0.9)
    ag.text(0.5 * sum(XH), 0.5, 'hole\nexit', color='0.3', ha='center', va='center')
    ag.set_xlim(*XL)
    ag.set_ylim(-0.02, 1.02)
    ag.set_xlabel(r'$x/D$', labelpad=1)
    ag.set_ylabel(r'$\gamma$ (mid-pitch)', labelpad=1)
    ag.text(0.012, 0.9, '(e)', transform=ag.transAxes, ha='left', va='center')
    ag.annotate('bands: $\\gamma$ with the detector\nthreshold halved (upper edge)\nand doubled (lower edge)',
                xy=band_xy, xytext=(2.7, 0.66), ha='left', va='center',
                arrowprops=dict(arrowstyle='->', color='0.3', lw=0.8))
    fig.savefig(OUT, dpi=600)
    print('wrote', OUT)


SPOTS_LAYOUT = os.environ.get('SPOTS_LAYOUT', 'reduced')       # 'reduced' (default): 193M, two plan views + (x, t) diagram; 'full': all panels


def main_reduced():
    """Reduced figure, 7.5 in x 3.1 in at scale 1: left, plan views of log10 d at the first and last snapshot time of the 193M (threshold
    contour, hole exit, mid-pitch line, tracked patch); right, (x, t) diagram on the mid-pitch line; one colour bar (top right)."""
    set_mpl_params(10.0 / 22.0, grid_on=False)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    HX, HY = load_hole()
    tag = 'c193b'
    lt = GV.policy_logth(tag)
    w0, sp, tk, Z, ld = pick_window(tag, lt)
    tf, frames = GV.plane_logd(tag, w0, w0 + TWIN)
    ik = [int(np.argmin(np.abs(tf - q))) for q in tk]
    tt = Z['tu']
    m = (tt >= w0) & (tt <= w0 + TWIN)
    line, tt = ld[m], tt[m] - w0
    shown = [np.asarray(frames[i], np.float32).ravel()[::7] for i in ik] + [ld[m].ravel()[::3]]
    vmin, vmax = np.percentile(np.concatenate(shown), [1, 99.5])        # percentiles of the 193M fields shown (the full figure pools 193M and LES)
    FH = 3.1
    fig = plt.figure(figsize=(FW, FH))
    f = lambda x0, y0, w, h: fig.add_axes([x0 / FW, y0 / FH, w / FW, h / FH])
    PW, PH, PG = 3.55, 1.18, 0.15
    xs = GV.XG
    im = None
    for row, k in enumerate((0, NSNAP - 1)):
        y0 = FH - 0.08 - (row + 1) * PH - row * PG
        ap = f(0.55, y0, PW, PH)
        yf, ff = fold_plane(GV.YG, np.asarray(frames[ik[k]], np.float32))
        im = ap.pcolormesh(xs, yf, ff.T, cmap='inferno', vmin=vmin, vmax=vmax, shading='nearest', rasterized=True)
        ap.contour(xs, yf, ff.T, levels=[lt], colors='w', linewidths=0.8, linestyles='solid', zorder=3)
        for dy in (0.0, 6.0):
            ap.add_patch(Polygon(np.c_[HX, HY + dy], closed=True, facecolor='w', edgecolor='k', lw=0.8, zorder=4))
        ap.plot([XL[0], XL[1]], [3.0, 3.0], color=OKABE_ITO[2], ls='--', lw=0.9, zorder=5)
        ap.set_xlim(*XL)
        ap.set_ylim(0, 6)
        ap.set_yticks([0, 3, 6])
        ap.set_ylabel(r'$y/D$', labelpad=1)
        ap.text(0.012, 0.9, r'$t_%d=%.1f$' % (row + 1, tf[ik[k]] - w0), transform=ap.transAxes, color='w', ha='left', va='top',
                zorder=6, bbox=BG)
        j = int(np.argmin(np.abs(sp['t'] - tf[ik[k]])))
        xc = 0.5 * (sp['le'][j] + sp['te'][j])
        ap.plot([xc], [2.35], marker='^', ms=5, mfc='w', mec='k', mew=0.7, ls='none', zorder=8)
        if row == 0:
            ap.set_xticklabels([])
        else:
            ap.set_xlabel(r'$x/D$', labelpad=1)
    XW, XH0, XB = 2.6, 1.9, 0.5
    ax = f(FW - XW - 0.16, XB, XW, XH0)
    xl = xs[:line.shape[1]]
    dtm = tt[1] - tt[0]
    ax.imshow(line, origin='lower', aspect='auto', cmap='inferno', vmin=vmin, vmax=vmax, interpolation='nearest',
              extent=[xl[0] - 0.01, xl[-1] + 0.01, tt[0] - dtm / 2, tt[-1] + dtm / 2])
    ax.contour(xl, tt, line, levels=[lt], colors='w', linewidths=0.8, linestyles='solid')
    for xh in XH:
        ax.axvline(xh, color='w', ls=':', lw=1.2)
    ax.set_xlim(*XL)
    ax.set_ylim(tt[0], tt[-1])
    ax.text(0.5 * sum(XH), tt[-1] - 0.03 * TWIN, 'hole\nexit', color='w', ha='center', va='top')
    tsn = tf[ik] - w0
    for r_, k in enumerate((0, NSNAP - 1)):                     # snapshot times: ticks on the left edge with the tags of the plan views
        ax.plot([XL[0], XL[0] + 0.7], [tsn[k], tsn[k]], color='w', lw=1.2, zorder=5)
        ax.text(XL[0] + 0.9, tsn[k], r'$t_%d$' % (r_ + 1), color='w', ha='left', va='top' if k == 0 else 'bottom',
                bbox=dict(facecolor='k', edgecolor='none', pad=0.5, alpha=0.85), zorder=6)
    ax.set_xlabel(r'$x/D$', labelpad=1)
    ax.set_ylabel(r'$(t-t_0)U_\infty/D$', labelpad=1)
    cax = f(FW - XW - 0.16, XB + XH0 + 0.3, XW, 0.1)
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    cax.xaxis.set_label_position('top')
    cb.set_label(r'$\log_{10}d$', labelpad=2)
    print(f'{tag}: threshold log10 d {lt:.3f}; window t U/D {w0:.1f}-{w0 + TWIN:.1f}; snapshot times (t-t0) {np.round(tsn, 1)}; '
          f'colour range {vmin:.2f} to {vmax:.2f}')
    fig.savefig(OUT, dpi=600)
    print('wrote', OUT, 'figsize in', FW, FH)


if __name__ == '__main__':
    main_reduced() if SPOTS_LAYOUT == 'reduced' else main()
