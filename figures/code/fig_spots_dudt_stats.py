#!/usr/bin/env python3
"""Helper of fig_spots (not a figure script): du/dt spot statistics on the z/D 0.1 mid-pitch data, x/D <= 20.

Method of the old fig_spots.py (git e0913c1): detector d = |du/dt| (u/U_inf, t in D/U) smoothed over 1.1 D/U, threshold 0.0385 U^2/D;
leading / trailing edge celerities from the x-t regions of the line y/D = 3 (edge_celerities, copied with the record-count rule
expressed in D/U); mergers and births from the band tracker (track_all, qualifies rules of asme_n2/n2_analysis.py) on the plan
detector cache, columns x/D <= 20 only.
A: old caches (detector from every record, 0.1 D/U) cut at x/D 20: line det_y+3.00 of <tag>_z010.npz, plane <tag>_z010_det2.npy.
B: the u lines cached for fig_spots_z010 (every second record, 0.2 D/U): detector recomputed, same threshold (line only;
   the plane tracker needs the plane, which is cached only in A).
usage: python3 fig_spots_dudt_stats.py   (output to stdout; log in the spots_z010 cache folder)
"""
import os
import sys
import json

import numpy as np
from scipy.ndimage import label, uniform_filter1d

SP = '/tmp/claude-1001/-home-sk999-Documents-ASME-2027/f4c92ceb-6df8-4254-96c3-6552d15d0ecf/scratchpad'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
sys.path.insert(0, SP)
import _gamma_extract as GX  # noqa: E402
import fig_spots_dudt_old as FS  # noqa: E402     # copy of git e0913c1:figures/code/fig_spots.py, kept beside this script

CACHE = '/mnt/Elements/DNS/analysis/asme2027_figs/spots_z010'
TH = 0.0385
U = 86.2812


def edge_celerities(line, t, x, minrec):
    m = line > TH
    lab, n = label(m, structure=np.ones((3, 3), int))
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
        use = (le - te) < 5.0
        bad = np.where(~use)[0]
        stop = bad[0] if len(bad) else len(rows)
        if stop < minrec:
            continue
        tt = t[rows[:stop]]
        if tt[-1] - tt[0] < 6.0:
            continue
        out.append((rows[0], t[rows[0]], te[0], np.polyfit(tt, le[:stop], 1)[0], np.polyfit(tt, te[:stop], 1)[0], tt[-1] - tt[0]))
    return np.array(out)


def stats(v):
    q = np.percentile(v, [25, 50, 75])
    return dict(median=round(float(q[1]), 3), q25=round(float(q[0]), 3), q75=round(float(q[2]), 3), mean=round(float(np.mean(v)), 3), n=int(len(v)))


def qualifies(L, life, t):
    if L['end'] is None or not (0.5 <= L['x0'] <= 8.0) or L['a0'] >= 0.1:
        return False
    tr = np.array(L['traj'])
    rend = L['end']
    pt = [q for q in life[L['into']]['traj'] if q[0] < rend]
    pa = pt[-1][2] if pt else 0
    if tr[:-1, 3].mean() > 0.02:
        return False
    life_d = t[int(tr[-1, 0])] - t[int(tr[0, 0])]
    return not (life_d < 6 or tr[:, 2].max() < 0.5 or pa < 0.5 or not (5 <= tr[-1, 1] <= 14))


def main():
    res = {}
    for tag in ('c193b', 'c20bw'):
        R = {}
        z = np.load(os.path.join(GX.CACHE, f'{tag}_z010.npz'))
        t = z['t']
        x = GX.XG
        ok = x <= 20.0 + 1e-9
        line = z['det_y+3.00'][:, ok]
        ec = edge_celerities(line, t, x[ok], 30)
        R['A_line_0.1'] = {nm: stats(ec[:, col]) for nm, col in (('leading', 3), ('trailing', 4))}
        R['A_line_0.1']['birth_x_median'] = float(np.median(ec[:, 2]))
        print(tag, 'A line', R['A_line_0.1'], flush=True)
        # B: lines of the 0.2 D/U records
        L = np.load(os.path.join(CACHE, f'{tag}_z010_lines.npz'))
        tb, xb, ub = L['t'], L['x'], L['u_p'].astype(np.float64) / U
        d = uniform_filter1d(np.abs(np.gradient(ub, tb, axis=0)), 5, axis=0, mode='nearest')
        okb = xb <= 20.0 + 1e-9
        ecb = edge_celerities(d[:, okb], tb, xb[okb], 15)
        R['B_line_0.2'] = {nm: stats(ecb[:, col]) for nm, col in (('leading', 3), ('trailing', 4))}
        R['B_line_0.2']['birth_x_median'] = float(np.median(ecb[:, 2]))
        print(tag, 'B line', R['B_line_0.2'], flush=True)
        # A plane tracker, x <= 20
        det = np.load(os.path.join(GX.CACHE, f'{tag}_z010_det2.npy'), mmap_mode='r')[:, :576, :]
        life = FS.track_all(det)
        xs, born = [], []
        for i, Lf in life.items():
            if qualifies(Lf, life, t):
                tr = np.array(Lf['traj'])
                xs.append(tr[-1, 1])
                born.append(Lf['x0'])
        R['A_mergers_x<=20'] = dict(stats(np.array(xs)), born_x_median=float(np.median(born))) if xs else {}
        R['identities'] = len(life)
        print(tag, 'A mergers', R['A_mergers_x<=20'], 'identities', len(life), flush=True)
        res[tag] = R
    json.dump(res, open(os.path.join(CACHE, 'dudt_stats.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
