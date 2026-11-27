#!/usr/bin/env python3
"""Shared extractor of the 3DNS / LES intermittency gamma, high-pass spanwise-stress detector (helper module, not a figure
script). Replaces the |du/dt| detector of _gamma_extract.py (kept unchanged for the record).

Detector (one definition for every quantity called intermittency in the paper):
  d = (v - <v>_T)^2 / U_inf^2, v the spanwise velocity (y direction), <v>_T the running mean in time over T = 2 D/U at every
  point (centred box, edge value repeated at the ends of the record), then d smoothed by a running mean of 1 D/U (same box).
  Laminar / turbulent: Otsu's method on log10 d (pooled histogram), gamma = fraction of time with d above the threshold.
  Indicator and Otsu classification after Dellacasagrande, Lengani, Bertolini, Sanz, J. Turbomach. 145(4) 041015 (2023), p. 3;
  the one change is the running mean in time in place of their spanwise mean (the flow is not homogeneous in span).
  No normalisation by gamma_max.
Time: the sampler records are unevenly spaced in the 193M run (0.08-0.38 D/U, crash repeat dropped by _gamma_extract.records).
Both runs are resampled linearly in time to the uniform step DT = 0.1 D/U (the median record spacing of both runs); the filter
windows are defined in D/U (T/DT, 1/DT records of the uniform series), not in sampler records.
Samplers: z010 plan view at z/D 0.1 (x/D -3..20.5 only, first NXR = 1176 columns of 1401, y/D -3..3 by 0.02, 301), x10 plane at
x/D 10 (y/D -3..3, 195 z points 0..3 D). Records rho u v w T float32 (v = column 2).
Windows: 193M C193b batch b010, 20M LES C20bW batch b006 (as _gamma_extract.py).
Cache /mnt/Elements/DNS/analysis/asme2027_figs/gamma_v/ (the staging copy of v and work files are on the local disk, STAGE):
  (STAGE) <tag>_<smp>_v.npy  staging copy of v on the sampler records (float32; deletable)
  <tag>_<smp>_v.npz        tu (uniform time, D/U), TS, LEV (log10 thresholds), hist (nT, NB) pooled log10 d histograms
                           (z010: x/D 0..20, all y; x10: all points), frac (nT, nLEV, ...) fraction of time with
                           log10 d > LEV for every point (z010: (x, y), x10: (y, z)), lines (z010, T = 2: log10 d at
                           y/D -3 and +3, (nt, NXR) each)
  <tag>_z010_logd.npy      log10 d (T = 2) of the whole plane, x/D -3..20.5, on every second uniform time step (0.2 D/U)
                           inside PLANE_WIN (60 D/U), float16; times in <tag>_z010_logd_t.npy
  frac is stored as uint16 (fraction * 65535); load() returns float32 fractions.
usage: python3 _gamma_extract_v.py c193b|c20bw z010|x10 [test]
"""
import json
import os
import sys
import time

import numpy as np
from scipy.ndimage import uniform_filter1d

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _gamma_extract as GX  # noqa: E402

CACHE = '/mnt/Elements/DNS/analysis/asme2027_figs/gamma_v'
# the staging copy of v (deletable, 5-7 GB) and the work files are written on the local disk (the NTFS volume of CACHE
# is slow under load); the finished result files are then copied to CACHE by publish().
STAGE = '/tmp/claude-1001/-home-sk999-Documents-ASME-2027/f4c92ceb-6df8-4254-96c3-6552d15d0ecf/scratchpad/gamma_v_stage'
PLANE_WIN = {'c193b': (140.0, 200.0), 'c20bw': (262.8, 322.8)}      # t U/D of the stored plan-view frames
D, U = GX.D, GX.U
RUNS = GX.RUNS
DT = 0.1
TS = (1.0, 2.0, 4.0)            # running-mean lengths T, D/U; TS[1] is the reference
TSM = 1.0                       # smoothing of d, D/U
NX, NXR, NY, NZ = GX.NX, 1176, GX.NY, GX.NZ
XG = GX.XG[:NXR]
YG = GX.YG
ZP = GX.ZP
NB, LO, HI = 1200, -12.0, 0.0
EDGES = np.linspace(LO, HI, NB + 1)
L0, LST, NL = -6.5, 0.05, 91
QS = 65535.0                    # frac is stored as uint16 = round(fraction * QS)
LEV = L0 + LST * np.arange(NL)
FLOOR = 1e-16


def logd_series(vu, nT, nsm):
    """vu (nt, ...) uniform-time v / U_inf -> log10 d."""
    vm = uniform_filter1d(vu, nT, axis=0, mode='nearest')
    d = uniform_filter1d((vu - vm) ** 2, nsm, axis=0, mode='nearest')
    return np.log10(np.maximum(d, FLOOR))


def resample_weights(tD):
    tu = tD[0] + DT * np.arange(int(np.floor((tD[-1] - tD[0]) / DT + 1e-9)) + 1)
    i1 = np.clip(np.searchsorted(tD, tu, 'right') - 1, 0, len(tD) - 2)
    w = np.clip((tu - tD[i1]) / (tD[i1 + 1] - tD[i1]), 0.0, 1.0)
    return tu, i1, w


def resample(v, i1, w):
    sh = (-1,) + (1,) * (v.ndim - 1)
    return (1 - w.reshape(sh)) * v[i1] + w.reshape(sh) * v[i1 + 1]


def exceed_fraction(ld, npts_shape):
    """ld (nt, ...) log10 d -> (NL, ...) fraction of time with ld above LEV[k]."""
    nt = ld.shape[0]
    flat = ld.reshape(nt, -1)
    b = np.clip(np.floor((flat - L0) / LST).astype(np.int64) + 1, 0, NL)       # d > LEV[k]  <=>  b >= k + 1
    npt = flat.shape[1]
    cnt = np.bincount((b + (NL + 1) * np.arange(npt)[None, :]).ravel(), minlength=npt * (NL + 1)).reshape(npt, NL + 1)
    ge = np.cumsum(cnt[:, ::-1], axis=1)[:, ::-1]                              # ge[:, j] = count with b >= j
    return np.rint(ge[:, 1:] / nt * QS).T.reshape((NL,) + tuple(npts_shape)).astype(np.uint16)


def hist_log(ld):
    i = np.clip(((ld - LO) / (HI - LO) * NB).astype(np.int64), 0, NB - 1)
    return np.bincount(i.ravel(), minlength=NB)


def stage_v(tag, smp, test=False):
    run, _ = RUNS[tag]
    os.makedirs(STAGE, exist_ok=True)
    out = os.path.join(STAGE, f'{tag}_{smp}_v.npy')
    idx, tD, info = GX.records(tag, smp)
    if test:
        idx, tD = idx[:3], tD[:3]
    json.dump(info, open(os.path.join(STAGE, f'{tag}_{smp}_info.json'), 'w'))
    np.save(os.path.join(STAGE, f'{tag}_{smp}_t.npy'), tD)
    done = out + '.done'
    if os.path.exists(done) and not test:
        return tD
    z = smp == 'z010'
    npts_f = NX * NY if z else NY * NZ
    npts_r = NXR * NY if z else NY * NZ
    shape = (len(idx), NXR, NY) if z else (len(idx), NY, NZ)
    mm = np.lib.format.open_memmap(out + ('.test' if test else ''), 'w+', np.float32, shape)
    fn = os.path.join(run, f'samp_{smp}.bin')
    assert idx.max() < os.path.getsize(fn) // (npts_f * 20)
    t0 = time.time()
    with open(fn, 'rb') as fh:
        for n, r in enumerate(idx):
            fh.seek(int(r) * npts_f * 20)
            q = np.fromfile(fh, '<f4', npts_r * 5).reshape(npts_r, 5)
            mm[n] = q[:, 2].reshape(shape[1:])
            if n % 100 == 0:
                print(f'{tag} {smp} read v record {n}/{len(idx)} [{time.time() - t0:.0f} s]', flush=True)
    mm.flush()
    if not test:
        open(done, 'w').write('ok')
    return tD


def detect(tag, smp, test=False):
    tD = stage_v(tag, smp, test)
    v = np.load(os.path.join(STAGE, f'{tag}_{smp}_v.npy' + ('.test' if test else '')), mmap_mode='r')
    tu, i1, w = resample_weights(tD)
    nt = len(tu)
    nT = [max(int(round(T / DT)), 1) for T in TS]
    nsm = max(int(round(TSM / DT)), 1)
    z = smp == 'z010'
    n0 = v.shape[1]
    step = 48 if z else n0
    frac = np.zeros((len(TS), NL) + v.shape[1:], np.uint16)
    hist = np.zeros((len(TS), NB))
    if z:
        sel2 = np.where((tu >= PLANE_WIN[tag][0]) & (tu <= PLANE_WIN[tag][1]))[0][::2]       # 0.2 D/U frames
        if test:
            sel2 = np.arange(0, nt, 2)
        logd2 = np.lib.format.open_memmap(os.path.join(STAGE, f'{tag}_z010_logd.npy' + ('.test' if test else '')), 'w+',
                                          np.float16, (len(sel2), NXR, NY))
        np.save(os.path.join(STAGE, f'{tag}_z010_logd_t.npy' + ('.test' if test else '')), tu[sel2])
        lines = {'y-3.00': np.zeros((nt, NXR), np.float32), 'y+3.00': np.zeros((nt, NXR), np.float32)}
    t0 = time.time()
    for i0 in range(0, n0, step):
        sl = slice(i0, min(i0 + step, n0))
        vu = resample(np.asarray(v[:, sl], np.float64), i1, w) / U
        for q, nTq in enumerate(nT):
            ld = logd_series(vu, nTq, nsm)
            frac[q, :, sl] = exceed_fraction(ld, ld.shape[1:])
            if z:
                ok = (XG[sl] >= -1e-9) & (XG[sl] <= 20.0 + 1e-9)
                if ok.any():
                    hist[q] += hist_log(ld[:, ok])
                if TS[q] == 2.0:
                    logd2[:, sl] = ld[sel2].astype(np.float16)
                    lines['y-3.00'][:, sl] = ld[:, :, 0]
                    lines['y+3.00'][:, sl] = ld[:, :, NY - 1]
            else:
                hist[q] += hist_log(ld)
        print(f'{tag} {smp} columns {i0}-{sl.stop}/{n0} [{time.time() - t0:.0f} s]', flush=True)
    if z:
        logd2.flush()
    extra = {f'line_{k}': a for k, a in lines.items()} if z else {}
    np.savez(os.path.join(STAGE, f'{tag}_{smp}_v{".test" if test else ""}.npz'), tu=tu, TS=np.array(TS), LEV=LEV, hist=hist,
             frac=frac, edges=EDGES, **extra)
    print(tag, smp, 'done: uniform records', nt, 'dt', DT, flush=True)
    if not test:
        publish(tag, smp)


def publish(tag, smp):
    """copy the result files of (tag, smp) from STAGE to CACHE (the staging copy of v stays local)."""
    import shutil
    os.makedirs(CACHE, exist_ok=True)
    names = [f'{tag}_{smp}_v.npz', f'{tag}_{smp}_t.npy', f'{tag}_{smp}_info.json']
    if smp == 'z010':
        names += [f'{tag}_z010_logd.npy', f'{tag}_z010_logd_t.npy']
    for n in names:
        t0 = time.time()
        shutil.copyfile(os.path.join(STAGE, n), os.path.join(CACHE, n + '.part'))
        os.replace(os.path.join(CACHE, n + '.part'), os.path.join(CACHE, n))
        print('copied', n, f'{time.time() - t0:.0f} s', flush=True)


def otsu(h, edges=EDGES):
    return GX.otsu(h, edges)


def load(tag, smp='z010'):
    f = os.path.join(STAGE, f'{tag}_{smp}_v.npz')          # the local copy first (fast); CACHE is the permanent copy
    if not os.path.exists(f):
        f = os.path.join(CACHE, f'{tag}_{smp}_v.npz')
    if not os.path.exists(f):
        raise SystemExit(f'{f} missing: run python3 _gamma_extract_v.py {tag} {smp}')
    Z = dict(np.load(f))
    Z['frac'] = Z['frac'].astype(np.float32) / QS
    return Z


def plane_logd(tag, t0, t1, T=2.0, pad=6.0):
    """log10 d of the whole z010 plane (x/D -3..20.5) for t0 <= t U/D <= t1 on every second uniform step (0.2 D/U), computed
    from the staging copy of v (the records within `pad` D/U of the window are filtered; the filters reach 1.1 D/U, so the
    window itself is identical to the full-record result). Returns (t, logd (nt, NXR, NY) float16)."""
    tD = np.load(os.path.join(STAGE, f'{tag}_z010_t.npy'))
    v = np.load(os.path.join(STAGE, f'{tag}_z010_v.npy'), mmap_mode='r')
    tu, i1, w = resample_weights(tD)
    sel = np.where((tu >= t0 - pad) & (tu <= t1 + pad))[0]
    lo, hi = i1[sel].min(), i1[sel].max() + 2
    vu = resample(np.asarray(v[lo:hi], np.float64), i1[sel] - lo, w[sel]) / U
    ld = logd_series(vu, max(int(round(T / DT)), 1), max(int(round(TSM / DT)), 1))
    keep = np.where((tu[sel] >= t0) & (tu[sel] <= t1))[0][::2]
    return tu[sel][keep], ld[keep].astype(np.float16)


def load_plane(tag):
    """log10 d (T = 2) of the plan view on 0.2 D/U frames inside PLANE_WIN: (t, logd (nt, NXR, NY) float16)."""
    out = []
    for n in (f'{tag}_z010_logd_t.npy', f'{tag}_z010_logd.npy'):
        f = os.path.join(STAGE, n)
        out.append(np.load(f if os.path.exists(f) else os.path.join(CACHE, n), mmap_mode='r'))
    return np.asarray(out[0]), out[1]


def gamma_at(frac, logth):
    """frac (NL, ...) -> gamma at the threshold 10**logth, linear interpolation in log10 threshold between LEV points."""
    k = (logth - L0) / LST
    k0 = int(np.clip(np.floor(k), 0, NL - 2))
    a = float(np.clip(k - k0, 0.0, 1.0))
    return (1 - a) * frac[k0] + a * frac[k0 + 1]


def line_logd(tag, T=2.0):
    """log10 d on the mid-pitch line y/D = +3 (sampler j = 300; j = 0 is identical), (nt, NXR), uniform step DT, from the
    staging copy of v."""
    tD = np.load(os.path.join(STAGE, f'{tag}_z010_t.npy'))
    v = np.load(os.path.join(STAGE, f'{tag}_z010_v.npy'), mmap_mode='r')
    tu, i1, w = resample_weights(tD)
    vu = resample(np.asarray(v[:, :, NY - 1], np.float64), i1, w) / U
    return logd_series(vu, max(int(round(T / DT)), 1), max(int(round(TSM / DT)), 1))


def build_thresholds():
    """Otsu thresholds (log10 d) for every T in TS, run, and pooling domain; written to thresholds.json (STAGE and CACHE).
    plane: z010 plane x/D 0..20, all y;  line: mid-pitch line y/D = 3, x/D 0..20."""
    import json as _json
    res = {}
    ok = (XG >= -1e-9) & (XG <= 20.0 + 1e-9)
    for q, T in enumerate(TS):
        r = {}
        for tag in ('c193b', 'c20bw'):
            Z = load(tag, 'z010')
            r[f'plane_{tag}'] = otsu(Z['hist'][q])
            ld = Z['line_y+3.00'] if T == 2.0 else line_logd(tag, T)
            h = hist_log(ld[:, ok])
            r[f'line_{tag}'] = dict(otsu(h), hist=h.tolist())
        res[f'{T:g}'] = r
    _json.dump(res, open(os.path.join(STAGE, 'thresholds.json'), 'w'))
    try:
        _json.dump(res, open(os.path.join(CACHE, 'thresholds.json'), 'w'))
    except OSError:
        pass
    return res


BAND_YLO = (1.5, 2.0, 2.5)      # band limits |y/D| >= ylo (to 3); 2.0 is the policy band


def band_hists(tag):
    """pooled histograms of log10 d over the bands |y/D| in [ylo, 3] (both sides of the pitch), x/D 0..20, for every T in TS and
    ylo in BAND_YLO -> {(T, ylo): hist}, from the staging copy of v."""
    tD = np.load(os.path.join(STAGE, f'{tag}_z010_t.npy'))
    v = np.load(os.path.join(STAGE, f'{tag}_z010_v.npy'), mmap_mode='r')
    tu, i1, w = resample_weights(tD)
    ix = np.where((XG >= -1e-9) & (XG <= 20.0 + 1e-9))[0]
    nsm = max(int(round(TSM / DT)), 1)
    H = {(T, yl): np.zeros(NB) for T in TS for yl in BAND_YLO}
    jy = {yl: np.abs(YG) >= yl - 1e-9 for yl in BAND_YLO}
    t0 = time.time()
    for c0 in range(0, len(ix), 100):
        sl = slice(ix[c0], ix[min(c0 + 100, len(ix)) - 1] + 1)
        vu = resample(np.asarray(v[:, sl], np.float64), i1, w) / U
        for T in TS:
            ld = logd_series(vu, max(int(round(T / DT)), 1), nsm)
            for yl in BAND_YLO:
                H[(T, yl)] += hist_log(ld[:, :, jy[yl]])
        print(f'{tag} band histograms columns {c0}/{len(ix)} [{time.time() - t0:.0f} s]', flush=True)
    return H


def add_band(res):
    import json as _json
    for tag in ('c193b', 'c20bw'):
        H = band_hists(tag)
        for (T, yl), h in H.items():
            res[f'{T:g}'][f'band{yl:g}_{tag}'] = dict(otsu(h), hist=h.tolist())
    for d in (STAGE, CACHE):
        try:
            _json.dump(res, open(os.path.join(d, 'thresholds.json'), 'w'))
        except OSError:
            pass
    return res


def thresholds():
    import json as _json
    f = os.path.join(STAGE, 'thresholds.json')
    if not os.path.exists(f):
        f = os.path.join(CACHE, 'thresholds.json')
    res = _json.load(open(f)) if os.path.exists(f) else build_thresholds()
    if 'band2_c193b' not in res['2']:
        res = add_band(res)
    return res


# Threshold policy (log10 d threshold of the Otsu classification):
#   'i'   each run its own Otsu threshold, pooled over the z010 plane x/D 0..20 (all y)
#   'ii'  the plane-pooled Otsu threshold of the 193M applied to both runs
#   'iii' the Otsu threshold of the 193M pooled over the mid-pitch line y/D = 3, x/D 0..20, applied to both runs
#   'band' (default) the Otsu threshold of the 193M pooled over the intermittent band |y/D| 2..3 (both sides of the pitch),
#          x/D 0..20, applied to both runs (owner decision 27 Nov)
#   'band1.5', 'band2.5'  the same with the band limits |y/D| 1.5..3 and 2.5..3 (sensitivity)
POLICY = 'band'


def policy_logth(tag, T=2.0, policy=None):
    """log10 of the threshold used for `tag` (see POLICY)."""
    p = policy or POLICY
    r = thresholds()[f'{T:g}']
    if p == 'i':
        return float(r[f'plane_{tag}']['log10'])
    if p == 'ii':
        return float(r['plane_c193b']['log10'])
    if p == 'iii':
        return float(r['line_c193b']['log10'])
    key = 'band2' if p == 'band' else p          # 'band1.5', 'band2.5'
    return float(r[key + '_c193b']['log10'])


FAC = (0.5, 2.0)        # threshold factors drawn as the shaded band around gamma(x) lines


def load_z010(tag, T=2.0, logth=None):
    """dict: gamma (NXR, NY) at 10**logth (default: policy_logth), x, y, tu, logth."""
    Z = load(tag, 'z010')
    q = int(np.argmin(np.abs(Z['TS'] - T)))
    if logth is None:
        logth = policy_logth(tag, T)
    return dict(gamma=gamma_at(Z['frac'][q], logth), x=XG, y=YG, tu=Z['tu'], logth=logth)


def load_x10(tag, T=2.0, logth=None):
    """dict: gamma (NY, NZ) at 10**logth, y, z."""
    Z = load(tag, 'x10')
    q = int(np.argmin(np.abs(Z['TS'] - T)))
    if logth is None:
        logth = policy_logth(tag, T)
    return dict(gamma=gamma_at(Z['frac'][q], logth), y=YG, z=ZP, logth=logth)


if __name__ == '__main__':
    detect(sys.argv[1], sys.argv[2], len(sys.argv) > 3)
