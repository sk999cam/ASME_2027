#!/usr/bin/env python3
"""Shared extractor of the 3DNS / LES intermittency gamma for fig_wall_state, fig_bl_profiles, fig_spots and fig_eta_x
(helper module, not a figure script).

Detector of analysis/layer_state (layer_samp.py, layer_whole_dns.py, FINDINGS LS1-LS6), copied unchanged:
  d = |du/dt|, u/U_inf and t in D/U, central difference on the record times (np.gradient), smoothed by a running mean of
  11 records (1.1 D/U, mode 'nearest'); gamma = fraction of records with d > threshold.
Thresholds: 0.0385 U^2/D (the one global Otsu threshold of the 50M, layer_otsu.py), 0.03 and 0.06 (sensitivity band).
Own Otsu: Otsu on the pooled histogram of log10 d (1200 bins, -6..1, d clipped at 1e-6; z010: x/D -3..20.5, all y; x10:
all y, z/D < 1.5), as layer_otsu.py / layer_whole_dns.py.

Samplers (same layout in every run): z010 plan view at z/D 0.1, x/D -3..25 by 0.02 (1401), y/D -3..3 by 0.02 (301), records
rho u v w T float32, point order x-major; x10 plane at x/D 10, y/D -3..3 by 0.02 (301) by 195 z points (0..3 D), same
record layout. Times: samp_z010.time column 3, in s; stat3_batches.txt (batch, duration, end time, ...) in s.
  193M C193b  batch b010 window = end - duration to end (4.413e-3 to 4.0342e-2 s). The sampler time file has a crash
              repeat (record 102, t goes back); the earlier copies of repeated times are dropped.
  20M LES     C20bW batch b006 window (6.560e-3 to 2.8117e-2 s).
Cache /mnt/Elements/DNS/analysis/asme2027_figs/gamma/: <tag>_z010_u.npy (staging copy of u, float32; deletable),
<tag>_z010.npz, <tag>_z010_det2.npy (detector on every second point, float16), <tag>_x10.npz.
usage: python3 _gamma_extract.py c193b|c20bw z010|x10      # staging + detection (u is read once from the run folder)
"""
import json
import os
import sys
import time

import numpy as np
from scipy.ndimage import uniform_filter1d

CACHE = '/mnt/Elements/DNS/analysis/asme2027_figs/gamma'
D, U = 0.00775, 86.2812
RUNS = {'c193b': ('/mnt/Elements/DNS/film_cooling/C193b_BR15', 10),
        'c20bw': ('/mnt/Tank/DNS/film_cooling/C20bW_BR15', 6)}
TWIN = 1.1
THS = (0.0385, 0.03, 0.06)
NX, NY, NZ = 1401, 301, 195
XG = -3.0 + 0.02 * np.arange(NX)
YG = -3.0 + 0.02 * np.arange(NY)
NB, LO, HI = 1200, -6.0, 1.0
EDGES = np.linspace(LO, HI, NB + 1)
LINE_J = {'y-3.00': 0, 'y-2.92': 4, 'y+2.92': 296, 'y+3.00': 300}
MAXRAM_GB = 20


def zprof(n_fine=60, z_fine=0.3, dz=0.02, z_top=3.0):
    zf = z_fine * (np.geomspace(1, 51, n_fine) - 1) / 50
    zc = np.arange(z_fine + dz, z_top + 1e-9, dz)
    return np.concatenate([zf, zc])


ZP = zprof()
assert len(ZP) == NZ


def window(tag):
    run, b = RUNS[tag]
    for ln in open(os.path.join(run, 'stat3_batches.txt')):
        f = ln.split()
        if int(f[0]) == b:
            return float(f[2]) - float(f[1]), float(f[2])
    raise SystemExit('no batch')


def records(tag, sampler='z010'):
    """indices into the sampler file and times (D/U) of the records inside the batch window, repeats removed."""
    run, _ = RUNS[tag]
    tt = np.loadtxt(os.path.join(run, f'samp_{sampler}.time'), ndmin=2)[:, 2]
    t0, t1 = window(tag)
    rmin = np.minimum.accumulate(tt[::-1])[::-1]
    later_min = np.concatenate([rmin[1:], [np.inf]])
    keep = tt < later_min                                  # drop a record if a later record has t <= its t
    inw = (tt >= t0 * (1 - 1e-12)) & (tt <= t1 * (1 + 1e-12))
    idx = np.where(keep & inw)[0]
    return idx, tt[idx] * U / D, dict(t0=t0, t1=t1, n_window=int(inw.sum()), n_dropped=int((~keep & inw).sum()),
                                      n_file=len(tt))


def detector(s, t, nwin):
    return uniform_filter1d(np.abs(np.gradient(s, t, axis=0)), nwin, axis=0, mode='nearest')


def hist_log(d):
    i = np.clip(((np.log10(np.maximum(d, 1e-6)) - LO) / (HI - LO) * NB).astype(np.int64), 0, NB - 1)
    return np.bincount(i.ravel(), minlength=NB)


def otsu(h, edges=EDGES):
    c = 0.5 * (edges[1:] + edges[:-1])
    p = h / h.sum()
    w0, m0 = np.cumsum(p), np.cumsum(p * c)
    mt = m0[-1]
    with np.errstate(divide='ignore', invalid='ignore'):
        sb = (mt * w0 - m0) ** 2 / (w0 * (1 - w0))
    sb[~np.isfinite(sb)] = -1
    k = int(np.argmax(sb))
    st2 = (p * (c - mt) ** 2).sum()
    return dict(log10=float(edges[k + 1]), th=float(10 ** edges[k + 1]), eta=float(sb[k] / st2), w0=float(w0[k]))


def stage_u(tag, sampler):
    """one pass over the run folder: u of every record in the window -> (nrec, ...) float32 npy."""
    run, _ = RUNS[tag]
    os.makedirs(CACHE, exist_ok=True)
    out = os.path.join(CACHE, f'{tag}_{sampler}_u.npy')
    idx, tD, info = records(tag, sampler)
    np.save(os.path.join(CACHE, f'{tag}_{sampler}_t.npy'), tD)
    json.dump(info, open(os.path.join(CACHE, f'{tag}_{sampler}_info.json'), 'w'))
    done = out + '.done'
    if os.path.exists(done):
        return idx, tD
    npts = NX * NY if sampler == 'z010' else NY * NZ
    shape = (len(idx), NX, NY) if sampler == 'z010' else (len(idx), NY, NZ)
    mm = np.lib.format.open_memmap(out, 'w+', np.float32, shape)
    fn = os.path.join(run, f'samp_{sampler}.bin')
    assert idx.max() < os.path.getsize(fn) // (npts * 20)
    t0 = time.time()
    with open(fn, 'rb') as fh:
        for n, r in enumerate(idx):
            fh.seek(int(r) * npts * 20)
            q = np.fromfile(fh, '<f4', npts * 5).reshape(npts, 5)
            mm[n] = q[:, 1].reshape(shape[1:])
            if n % 200 == 0:
                print(f'{tag} {sampler} read record {n}/{len(idx)} [{time.time() - t0:.0f} s]', flush=True)
    mm.flush()
    open(done, 'w').write('ok')
    return idx, tD


def detect_z010(tag):
    idx, tD = stage_u(tag, 'z010')
    u = np.load(os.path.join(CACHE, f'{tag}_z010_u.npy'), mmap_mode='r')
    nr = len(tD)
    nwin = int(round(TWIN / np.median(np.diff(tD))))
    gam = np.zeros((len(THS), NX, NY), np.float32)
    um, ur = np.zeros((NX, NY)), np.zeros((NX, NY))
    h = np.zeros(NB)
    det2 = np.lib.format.open_memmap(os.path.join(CACHE, f'{tag}_z010_det2.npy'), 'w+', np.float16,
                                     (nr, (NX + 1) // 2, (NY + 1) // 2))
    lines = {k: np.zeros((nr, NX), np.float32) for k in LINE_J}
    t0 = time.time()
    for i0 in range(0, NX, 100):
        c = min(100, NX - i0)
        s = np.asarray(u[:, i0:i0 + c], np.float64) / U
        um[i0:i0 + c], ur[i0:i0 + c] = s.mean(0), s.std(0)
        d = detector(s, tD, nwin)
        del s
        for q, th in enumerate(THS):
            gam[q, i0:i0 + c] = (d > th).mean(0)
        ok = XG[i0:i0 + c] <= 20.5 + 1e-9
        h += hist_log(d[:, ok])
        for k, j in LINE_J.items():
            lines[k][:, i0:i0 + c] = d[:, :, j]
        sel = np.arange(i0, i0 + c)
        sel = sel[sel % 2 == 0]
        det2[:, sel // 2] = d[:, sel - i0, ::2].astype(np.float16)
        print(f'{tag} z010 x columns {i0}-{i0 + c} [{time.time() - t0:.0f} s]', flush=True)
    det2.flush()
    o600 = otsu(h.reshape(600, 2).sum(1), np.linspace(LO, HI, 601))
    o1200 = otsu(h)
    np.savez(os.path.join(CACHE, f'{tag}_z010.npz'), gamma=gam, TH=np.array(THS), t=tD, hist=h, umean=um, urms=ur,
             nwin=nwin, otsu_th=o1200['th'], otsu_th600=o600['th'], otsu_eta=o1200['eta'],
             **{f'det_{k}': v for k, v in lines.items()})
    print(tag, 'z010 records', nr, 'nwin', nwin, 'Otsu', o1200, 'Otsu600', o600)


def detect_x10(tag):
    idx, tD = stage_u(tag, 'x10')
    u = np.load(os.path.join(CACHE, f'{tag}_x10_u.npy'))
    nwin = int(round(TWIN / np.median(np.diff(tD))))
    d = detector(u.astype(np.float64) / U, tD, nwin)
    gam = np.stack([(d > th).mean(0) for th in THS]).astype(np.float32)
    zl = ZP < 1.5
    h_lt, h_all = hist_log(d[:, :, zl]), hist_log(d)
    o_lt, o_all = otsu(h_lt), otsu(h_all)
    np.savez(os.path.join(CACHE, f'{tag}_x10.npz'), gamma=gam, TH=np.array(THS), t=tD, z=ZP, y=YG, nwin=nwin,
             hist_z15=h_lt, hist_all=h_all, otsu_th_z15=o_lt['th'], otsu_th_all=o_all['th'])
    print(tag, 'x10 records', len(tD), 'Otsu z<1.5', o_lt, 'all z', o_all)


def load_z010(tag):
    f = os.path.join(CACHE, f'{tag}_z010.npz')
    if not os.path.exists(f):
        detect_z010(tag)
    return dict(np.load(f))


def load_x10(tag):
    f = os.path.join(CACHE, f'{tag}_x10.npz')
    if not os.path.exists(f):
        detect_x10(tag)
    return dict(np.load(f))


if __name__ == '__main__':
    tag, smp = sys.argv[1], sys.argv[2]
    (detect_z010 if smp == 'z010' else detect_x10)(tag)
