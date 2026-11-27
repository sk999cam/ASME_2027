#!/usr/bin/env python3
"""Helper of fig_spots_z010 (not a figure script): mid-pitch lines of the z/D 0.1 plane sampler.

For every record of the batch window whose time matches a record of the true-wall-C_f set (n2_spots_cf/<tag>_spacetime.npz,
record spacing 0.2 D/U, i.e. every second sampler record) the lines y/D = -3.00 (j = 0) and y/D = +3.00 (j = 300) are
kept for x/D -3 to 21 (ix 0..1200), fields u, T, rho (float32, SI). The sampler layout is that of _gamma_extract.py:
records of 1401 x 301 points (x-major, y fastest) with rho u v w T float32; only the first 1201 x-columns of every record
are read (7.2 MB of 8.4 MB). Records are selected with the window logic of _gamma_extract.records.
Cache: /mnt/Elements/DNS/analysis/asme2027_figs/spots_z010/<tag>_z010_lines.npz (+ .part progress files).
usage: python3 _spots_z010_extract.py c193b|c20bw [test]
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _gamma_extract as GX  # noqa: E402

OUT = '/mnt/Elements/DNS/analysis/asme2027_figs/spots_z010'
N2 = '/mnt/Elements/DNS/analysis/asme2027_figs/n2_spots_cf'
NXR = 1201


def select(tag):
    idx, tD, info = GX.records(tag, 'z010')
    tc = np.load(os.path.join(N2, f'{tag}_spacetime.npz'))['t']
    j = np.abs(tD[:, None] - tc[None, :]).argmin(0)
    err = np.abs(tD[j] - tc)
    assert err.max() < 0.08, err.max()      # the 193M sampler spacing is not uniform (0.077 to 0.199 D/U)
    return idx[j], tD[j], err.max()


def main(tag, test=False):
    os.makedirs(OUT, exist_ok=True)
    run, _ = GX.RUNS[tag]
    idx, tD, err = select(tag)
    if test:
        idx, tD = idx[:2], tD[:2]
    n = len(idx)
    npts = NXR * GX.NY
    print(f'{tag}: {n} records, max time mismatch to the C_f set {err:.4f} D/U, read {n * npts * 20 / 1e9:.2f} GB', flush=True)
    A = {f'{f}_{s}': np.zeros((n, NXR), np.float32) for f in ('u', 'T', 'rho') for s in ('m', 'p')}
    t0 = time.time()
    with open(os.path.join(run, 'samp_z010.bin'), 'rb') as fh:
        for k, r in enumerate(idx):
            fh.seek(int(r) * GX.NX * GX.NY * 20)
            q = np.fromfile(fh, '<f4', npts * 5).reshape(NXR, GX.NY, 5)
            for s, j in (('m', 0), ('p', GX.NY - 1)):
                A[f'rho_{s}'][k], A[f'u_{s}'][k], A[f'T_{s}'][k] = q[:, j, 0], q[:, j, 1], q[:, j, 4]
            if k % 100 == 0:
                print(f'{tag} record {k}/{n} [{time.time() - t0:.0f} s]', flush=True)
    if test:
        for key, v in A.items():
            print(key, v[0, :3], v[0, 200], v[0, 800])
        return
    np.savez(os.path.join(OUT, f'{tag}_z010_lines.npz'), t=tD, x=GX.XG[:NXR], rec=idx, **A)
    print(tag, 'done', f'{time.time() - t0:.0f} s', flush=True)


if __name__ == '__main__':
    main(sys.argv[1], len(sys.argv) > 2)
