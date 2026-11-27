"""Read-only check of the intermittency numbers (Section 4.1, layer between the jets). Run with python3 -B. Writes nothing."""
import sys

import numpy as np

sys.path.insert(0, '/home/sk999/Documents/ASME_2027/figures/code')
import _gamma_extract_v as GV  # noqa: E402

xg, yg = GV.XG, GV.YG
jm = np.where(np.abs(yg) >= 2.81 - 1e-9)[0]
print('policy', GV.POLICY, 'x range', xg.min(), xg.max(), 'y range', yg.min(), yg.max())
out = {}
for tag in ('c193b', 'c20bw'):
    G = GV.load_z010(tag)
    g = G['gamma']
    gm = np.nanmean(g[:, jm], axis=1)
    out[tag] = gm
    ok = xg <= 20 + 1e-9
    first = xg[ok][np.where(gm[ok] > 0.8)[0][0]] if np.any(gm[ok] > 0.8) else None
    sel = (xg >= 0) & (xg <= 20)
    print(f'{tag}: logth {G["logth"]:.3f}; gamma_mid first > 0.8 at x/D {first}; gamma_mid at x/D 20 (19.74-20.26) '
          f'{np.nanmean(gm[np.abs(xg - 20) < 0.26]):.3f}; at 19.5-20.0 {np.nanmean(gm[(xg >= 19.5) & (xg <= 20)]):.3f};'
          f' plane mean x/D 0-20 {np.nanmean(g[sel]):.3f}; x/D -3-20 {np.nanmean(g[xg <= 20]):.3f}')
d = out['c20bw'] - out['c193b']
for lo in (0, 5, 10, 15):
    s = (xg >= lo) & (xg <= 20)
    print(f'LES - 193M gamma_mid, x/D {lo}-20: max |diff| {np.nanmax(np.abs(d[s])):.3f}')
try:
    Z = np.load('/mnt/Elements/DNS/analysis/asme2027_figs/bl_extract/rans50_z010_cells.npz')
    print('RANS file keys', list(Z.keys()))
except Exception as e:  # noqa: BLE001
    print('RANS file', e)
