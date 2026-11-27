#!/usr/bin/env python3
"""Helper of fig_spots (not a figure script): spot celerities with the high-pass spanwise-stress detector.

Method of fig_spots_dudt_stats.py (old |du/dt| detector, kept for the record): on the binary mid-pitch (x, t) field
(line y/D = +3, log10 d above the threshold, x/D <= 20, uniform step 0.1 D/U) every 8-connected region is a candidate spot;
it is kept if it is born small (width at its first record <= 0.6 D, centre between x/D 0 and 14), is followed while its
width is below 5 D, lasts at least 6 D/U and at least 30 records (3 D/U); leading / trailing edge celerity = least-squares
slope of the largest / smallest x of the region against time. Prints median, quartiles, mean and N of both celerities per run
and threshold policy (iii: 193M mid-pitch-line Otsu threshold for both runs, the policy of the figures, see _gamma_extract_v.py;
ii: 193M plane-pooled threshold for both runs; i: own plane-pooled threshold).
usage: python3 fig_spots_stats_v.py
"""
import json
import os
import sys

import numpy as np
from scipy.ndimage import label

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _gamma_extract_v as GV  # noqa: E402

MINREC = 30


def find_spots(mask, t, x, minrec=MINREC):
    """mask (nt, nx) boolean on the line, t, x -> list of dicts (row0, t0, x0, cle, cte, life, t, le, te)."""
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
        use = (le - te) < 5.0
        bad = np.where(~use)[0]
        stop = bad[0] if len(bad) else len(rows)
        if stop < minrec:
            continue
        tt = t[rows[:stop]]
        if tt[-1] - tt[0] < 6.0:
            continue
        out.append(dict(row0=int(rows[0]), t0=float(t[rows[0]]), x0=float(te[0]), cle=float(np.polyfit(tt, le[:stop], 1)[0]),
                        cte=float(np.polyfit(tt, te[:stop], 1)[0]), life=float(tt[-1] - tt[0]), t=tt, le=le[:stop],
                        te=te[:stop]))
    return out


def stats(v):
    q = np.percentile(v, [25, 50, 75])
    return dict(median=round(float(q[1]), 3), q25=round(float(q[0]), 3), q75=round(float(q[2]), 3),
                mean=round(float(np.mean(v)), 3), n=int(len(v)))


def line(tag):
    Z = GV.load(tag, 'z010')
    ok = GV.XG <= 20.0 + 1e-9
    return Z['tu'], GV.XG[ok], Z['line_y+3.00'][:, ok]


def main():
    res = {}
    for pol in ('iii', 'ii', 'i'):
        GV.POLICY = pol
        for tag in ('c193b', 'c20bw'):
            t, x, ld = line(tag)
            lt = GV.policy_logth(tag)
            sp = find_spots(ld > lt, t, x)
            r = {'log10_th': round(float(lt), 3), 'N': len(sp)}
            if sp:
                r['leading'] = stats(np.array([s['cle'] for s in sp]))
                r['trailing'] = stats(np.array([s['cte'] for s in sp]))
                r['birth_x_median'] = round(float(np.median([s['x0'] for s in sp])), 2)
                r['life_median_DU'] = round(float(np.median([s['life'] for s in sp])), 1)
            res[f'{tag}_policy_{pol}'] = r
            print(tag, 'policy', pol, json.dumps(r), flush=True)
    json.dump(res, open(os.path.join(GV.CACHE if os.path.isdir(GV.CACHE) else GV.STAGE, 'spot_celerities_v.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
