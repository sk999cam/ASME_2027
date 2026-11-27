#!/usr/bin/env python3
"""F9 fig_spots: growth and merger of turbulent spots between the jets, 193M 3DNS, plane sampler z010 (z/D 0.1), b010 window.

(a) Space-time diagram of the detector d = |du/dt| (smoothed over 1.1 D/U, as the intermittency gamma of F8) along the
    mid-pitch sampler line y/D = +3 (the periodic line, same as y/D = -3): x/D -3 to 20.5 against t U_inf/D; white contour at
    the global Otsu threshold 0.0385 U^2/D; the dashed vertical lines mark the times of the snapshots in (b); the light-blue
    curve is the streamwise centroid of the spot tracked in (b).
(b) Six plan views of d at chosen times (y/D 0 to 6, mid-pitch y/D = 3 in the centre, the jet wakes at the edges, folded through
    the periodic boundary), with the Otsu contour (white) and the outline of ONE spot, tracked record by record from its birth
    to its merger with a neighbour and after it (light blue).
Spot tracking: in each record the patches of d > 0.0385 U^2/D in the between-jet band |y/D| >= 2 (8-connected, area >= 0.02 D^2,
    as analysis/layer_state/layer_samp.py --spots) are linked to the patches of the previous record by overlap; the patch
    that overlaps the largest earlier patch keeps its identity, the others end there (merger). A spot is the identity born as a
    patch with no earlier overlap. The shown spot is chosen automatically: born at x/D 0.5 to 8 as a patch smaller than 0.1 D^2,
    at least 6 D/U old, 0.5 D^2 large and clear of the band edge |y/D| = 2 (touching it in at most 2 % of the records) before it ends in a merger with a patch of at least 0.5 D^2 at x/D 5 to 14.
Celerities (printed): in the x-t diagram of the line y/D = 3 each connected region (8-connected, d > 0.0385) that starts
    narrower than 0.6 D at x/D 0 to 14 is one spot cut by the line; while it is narrower than 5 D, its leading edge (largest x)
    and trailing edge (smallest x) are fitted by least squares x(t); the slope is the celerity over U_inf. Median and
    inter-quartile range over the spots with at least 6 D/U of fit.
Also printed: gamma_mid(x), the time fraction of d > threshold averaged over the band |y/D| 2.81 to 3.0 and at y/D 3 alone; the
    spot birth rate from the tracker (identities born at x/D 0 to 6 in the band that live at least 5 D/U).
Data: /mnt/Elements/DNS/film_cooling/C193b_BR15 samp_z010, the records inside the stat3 batch b010 window, via _gamma_extract.py
    (cache /mnt/Elements/DNS/analysis/asme2027_figs/gamma/).
usage: python3 fig_spots.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from cycler import cycler  # noqa: E402
from scipy.ndimage import label  # noqa: E402

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from functions import set_mpl_params  # noqa: E402
import _gamma_extract as GX  # noqa: E402

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_spots.png'
TAG = 'c193b'
TH = 0.0385
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
MARK = '#56B4E9'
XG2 = -3.0 + 0.04 * np.arange(701)
YG2 = -3.0 + 0.04 * np.arange(151)
DA = 0.04 ** 2
AMIN = 0.02
JA = np.where(YG2 >= 2 - 1e-9)[0]                       # y 2 .. 3
JB = np.where(YG2 <= -2 + 1e-9)[0][1:]                  # y -2.96 .. -2 -> y + 6 = 3.04 .. 4
YB = np.concatenate([YG2[JA], YG2[JB] + 6.0])


def band_labels(f):
    m = np.concatenate([f[:, JA], f[:, JB]], 1) > TH
    lab, n = label(m, structure=np.ones((3, 3), int))
    if n:
        sz = np.bincount(lab.ravel())
        small = np.where(sz * DA < AMIN)[0]
        small = small[small > 0]
        if len(small):
            lab[np.isin(lab, small)] = 0
    return lab


def track_all(det):
    """identity of every patch through all records -> dict id: dict(r0, x0, a0, traj[(r, xc, area)], end, into, merges)."""
    life = {}
    prev = None
    nid = 0
    for r in range(det.shape[0]):
        lab = band_labels(np.asarray(det[r], np.float32))
        ids = [k for k in np.unique(lab) if k > 0]
        pids = np.zeros(lab.max() + 1, int)
        for k in ids:
            mk = lab == k
            area = mk.sum() * DA
            xc = XG2[np.nonzero(mk)[0]].mean()
            edge = bool(mk[:, 0].any() or mk[:, -1].any())
            pred = []
            if prev is not None:
                pl, pi = prev
                pred = sorted({pi[o] for o in np.unique(pl[mk]) if o > 0})
            if not pred:
                nid += 1
                my = nid
                life[my] = dict(r0=r, x0=xc, a0=area, traj=[], end=None, into=None, merges=[])
            else:
                my = max(pred, key=lambda p: life[p]['traj'][-1][2])
                for p in pred:
                    if p != my and life[p]['end'] is None:
                        life[p]['end'], life[p]['into'] = r, my
                        life[my]['merges'].append((r, xc, life[p]['traj'][-1][2], area))
            life[my]['traj'].append((r, xc, area, edge))
            pids[k] = my
        prev = (lab, pids)
    return life


def choose(life, t):
    best = None
    for i, L in life.items():
        if L['end'] is None or not (0.5 <= L['x0'] <= 8.0) or L['a0'] >= 0.1:
            continue
        tr = np.array(L['traj'])
        rend = L['end']
        life_d = t[int(tr[-1, 0])] - t[int(tr[0, 0])]
        # partner: the patch that keeps its identity at the merger, its area in the record before
        pt = [q for q in life[L['into']]['traj'] if q[0] < rend]
        pa = pt[-1][2] if pt else 0
        if tr[:-1, 3].mean() > 0.02:        # the patch touches the band edge |y/D| = 2: not a free spot
            continue
        if life_d < 6 or tr[:, 2].max() < 0.5 or pa < 0.5 or not (5 <= tr[-1, 1] <= 14):
            continue
        score = life_d
        if best is None or score > best[0]:
            best = (score, i)
    return best[1] if best else None


def follow(det, life, i, r_end):
    """per-record masks (full 701 x 151 plan) of spot i from its birth to r_end, continuing through the merger."""
    L = life[i]
    r0 = L['r0']
    masks = {}
    lab = band_labels(np.asarray(det[r0], np.float32))
    xs = np.array([[XG2[a], 0] for a in range(0)])
    best, bd = None, 1e9
    for k in np.unique(lab):
        if k == 0:
            continue
        mk = lab == k
        d = abs(XG2[np.nonzero(mk)[0]].mean() - L['x0'])
        if d < bd:
            best, bd = k, d
    cur = lab == best
    for r in range(r0, r_end + 1):
        lab = band_labels(np.asarray(det[r], np.float32))
        if r > r0:
            ov = np.bincount(lab[cur], minlength=lab.max() + 1)
            ov[0] = 0
            if ov.max() == 0:
                break
            cur = lab == int(np.argmax(ov))
        masks[r] = cur.copy()
    return masks


def to_plan(mask):
    """band mask (701 x len(YB)) -> plan array (701 x 151) in the folded coordinate y' = y (y >= 0), y + 6 (y < 0)."""
    full = np.zeros((701, 151), bool)
    full[:, JA] = mask[:, :len(JA)]
    full[:, JB] = mask[:, len(JA):]
    return full


def fold(a):
    """plan array on y -3..3 -> folded y' 0..6 (mid-pitch y' = 3 in the centre)."""
    pos = a[:, 75:]                      # y 0 .. 3
    neg = a[:, 1:75]                     # y -2.96 .. -0.04 -> y' 3.04 .. 5.96
    return np.concatenate([pos, neg], 1), np.concatenate([YG2[75:], YG2[1:75] + 6.0])


def edge_celerities(det_line, t, x):
    """leading and trailing edge slopes of the x-t regions on the line."""
    m = det_line > TH
    lab, n = label(m, structure=np.ones((3, 3), int))
    out = []
    for k in range(1, n + 1):
        rr, ii = np.nonzero(lab == k)
        rows = np.unique(rr)
        if len(rows) < 30:
            continue
        le = np.array([x[ii[rr == r]].max() for r in rows])
        te = np.array([x[ii[rr == r]].min() for r in rows])
        if (le[0] - te[0]) > 0.6 or not (0.0 <= 0.5 * (le[0] + te[0]) <= 14.0):
            continue
        use = (le - te) < 5.0
        # contiguous from the start
        bad = np.where(~use)[0]
        stop = bad[0] if len(bad) else len(rows)
        if stop < 30:
            continue
        tt = t[rows[:stop]]
        if tt[-1] - tt[0] < 6.0:
            continue
        cl = np.polyfit(tt, le[:stop], 1)[0]
        ct = np.polyfit(tt, te[:stop], 1)[0]
        out.append((rows[0], t[rows[0]], te[0], cl, ct, tt[-1] - tt[0]))
    return np.array(out)


def main():
    z = GX.load_z010(TAG)
    t = z['t']
    det = np.load(os.path.join(GX.CACHE, f'{TAG}_z010_det2.npy'), mmap_mode='r')
    x = GX.XG
    line = z['det_y+3.00']
    line_m = z['det_y-3.00']
    print('records', len(t), 't', t[0], t[-1], 'D/U')
    print('max |d(y=+3) - d(y=-3)|', np.abs(line - line_m).max(), '(periodic line check)')
    G = z['gamma'][0]
    okx = x <= 20.5 + 1e-9
    jm = np.where(np.abs(GX.YG) >= 2.81 - 1e-9)[0]
    print('gamma_mid (band |y| 2.81-3.0 | line y=+3), threshold 0.0385:')
    for q in (2, 4, 6, 8, 10, 12, 14, 16, 18, 20):
        k = int(round((q + 3) / 0.02))
        print(f'  x/D {q}: {G[k, jm].mean():.3f} | {G[k, 300]:.3f}')

    # --- celerities from (a)
    ec = edge_celerities(line[:, okx], t, x[okx])
    if len(ec):
        for nm, col in (('leading edge', 3), ('trailing edge', 4)):
            v = ec[:, col]
            q = np.percentile(v, [25, 50, 75])
            print(f'celerity {nm}: median {q[1]:.2f} U_inf, IQR {q[0]:.2f}-{q[2]:.2f}, mean {v.mean():.2f}, N {len(v)}')
        print('fit duration median', np.median(ec[:, 5]), 'D/U; birth x/D median', np.median(ec[:, 2]))

    # --- tracking
    life = track_all(det)
    sel = choose(life, t)
    print('tracked identities', len(life), 'chosen', sel)
    L = life[sel]
    tr = np.array(L['traj'])
    r_b, r_e = L['r0'], L['end']
    print(f'spot: born t {t[r_b]:.1f} D/U x/D {L["x0"]:.2f}; merges at t {t[r_e]:.1f} D/U x/D {tr[-1, 1]:.2f}; max area {tr[:, 2].max():.2f} D^2')
    nb = sum(1 for i, a in life.items() if 0 <= a['x0'] <= 6 and len(a['traj']) * 0.1 >= 5.0)
    area_band = 6.0 * 2.0 / 2 * 1.0   # band |y| >= 2 is 2 D wide over x 0..6: 12 D^2
    print(f'births x/D 0-6 (live >= 5 D/U): {nb} in {t[-1] - t[0]:.0f} D/U over 12 D^2 -> {nb / 12.0 / (t[-1] - t[0]) * 100:.3f} per D^2 per 100 D/U')
    r_after = min(r_e + int(round(8 / np.median(np.diff(t)))), len(t) - 1)
    masks = follow(det, life, sel, r_after)
    rs_all = sorted(masks)
    # snapshot records: birth, 3 growth, just before merger, just after
    r_pre, r_post = r_e - 2, min(r_e + 10, rs_all[-1])
    mid = [int(round(a)) for a in np.linspace(r_b, r_pre, 5)[1:4]]
    rsnap = [r_b + 3, mid[0], mid[1], mid[2], r_pre, r_post]
    rsnap = [r for r in rsnap if r in masks]
    print('snapshots (record, t U/D, x/D of the spot centroid):')
    for r in rsnap:
        mk = masks[r]
        xc = XG2[np.nonzero(mk)[0]].mean()
        print(f'  {r}  t {t[r]:.1f}  x {xc:.2f}  area {mk.sum() * DA:.2f} D^2')

    # --- figure
    set_mpl_params(0.7, grid_on=False)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    fig = plt.figure(figsize=(15.0, 10.5))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.9, 1, 1], hspace=0.28, wspace=0.12)
    axa = fig.add_subplot(gs[0, :])
    norm = LogNorm(1e-2, 3e-1)
    tt = np.concatenate([[t[0] - 0.05], 0.5 * (t[1:] + t[:-1]), [t[-1] + 0.05]])
    xe = np.concatenate([[x[0] - 0.01], 0.5 * (x[1:] + x[:-1]), [x[-1] + 0.01]])
    ke = np.searchsorted(xe, 20.5)
    im = axa.pcolormesh(tt, xe[:ke + 1], np.maximum(line[:, :ke].T, 1e-4), cmap='inferno', norm=norm, shading='flat', rasterized=True)
    axa.contour(t, x[:ke], line[:, :ke].T, levels=[TH], colors='w', linewidths=0.5)
    axa.plot(t[tr[:, 0].astype(int)], tr[:, 1], color=MARK, lw=1.6)
    cols = ['#56B4E9']
    for r in rsnap:
        axa.axvline(t[r], color=MARK, ls='--', lw=0.8)
    axa.set_xlim(t[0], t[-1])
    axa.set_ylim(-3, 20.5)
    axa.set_xlabel(r'$tU_\infty/D$')
    axa.set_ylabel(r'$x/D$')
    axa.text(0.006, 0.93, r'(a) $y/D=3$', transform=axa.transAxes, color='k', ha='left', va='center', bbox=dict(fc='w', ec='none', pad=1.5))
    cb = fig.colorbar(im, ax=axa, pad=0.01, fraction=0.02)
    cb.set_label(r'$|\partial u/\partial t|$ smoothed, $U_\infty^2/D$')
    cb.ax.axhline(TH, color='w', lw=1)

    for n, r in enumerate(rsnap):
        a = fig.add_subplot(gs[1 + n // 3, n % 3])
        f = np.asarray(det[r], np.float32)
        ff, yf = fold(f)
        a.pcolormesh(XG2, yf, np.maximum(ff.T, 1e-4), cmap='inferno', norm=norm, shading='auto', rasterized=True)
        a.contour(XG2, yf, ff.T, levels=[TH], colors='w', linewidths=0.4)
        mk, _ = fold(to_plan(masks[r]).astype(float))
        a.contour(XG2, yf, mk.T, levels=[0.5], colors=MARK, linewidths=1.6)
        a.set_xlim(-3, 20.5)
        a.set_ylim(0, 6)
        a.set_yticks([0, 3, 6])
        a.set_aspect('auto')
        a.text(0.012, 0.86, rf'({chr(98 + n)}) $tU_\infty/D={t[r]:.0f}$', transform=a.transAxes, color='k', ha='left', va='center', bbox=dict(fc='w', ec='none', pad=1.5))
        if n // 3 == 1:
            a.set_xlabel(r'$x/D$')
        if n % 3 == 0:
            a.set_ylabel(r'$y/D$')
    fig.savefig(OUT, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    main()
