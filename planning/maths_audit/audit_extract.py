#!/usr/bin/env python3
"""Extract the point set of the RANS-50M vs 193M-b010 comparison (as meff_rans50_c193b10) into one npz.

Read-only on all existing code and data. Imports errorbars.py (analysis/meff_c50b_c100b) for the point set, the
193M interpolation and the RANS reader, so the points are those of regions_se.txt. Writes only to the scratchpad.

Per point (3.14 M, every 2nd 50M grid point of the point-wise regions):
  reg, w (= 8 dV), 193M: rho k L cos S (interpolated, full record b010);
  RANS: rho k L mut (L = mut/(rho k^1/2)), S_r (own strain);
  50M 3DNS b006 (same grid, no interpolation): rho k L cos S meff, eps (W/m3, budgets_c50b/b006).
usage: PYTHONDONTWRITEBYTECODE=1 MEFF_REGIONS=pp MEFF_SFX_C193b=_b010 python3 -B audit_extract.py OUT.npz
"""
import os
import sys

import numpy as np

os.environ['MEFF_REGIONS'] = 'pp'
os.environ['MEFF_SFX_C193b'] = '_b010'
sys.dont_write_bytecode = True
EB = '/home/sk999/Documents/4dns_filmcooled/analysis/meff_c50b_c100b'
sys.path.insert(0, EB)
import errorbars as eb  # noqa: E402

out = sys.argv[1]
names, P = eb.points(eb.CASES['C50b'], 2)
print(len(P['n']), 'points', names, flush=True)
cache = os.path.join(eb.W, 'locate_C50b_points_in_C193b_pp.npz')
assert int(np.load(cache)['n']) == len(P['xyz']), 'locate cache does not match the point set; would write to Tank'
ref = eb.at_points_ref('C193b', P['xyz'], cache, False)
print('193M read', flush=True)
R = eb.rans_at(P, 'RANS')
N = len(P['n'])
d50 = {q: np.full(N, np.nan) for q in ('rho_d', 'k_d', 'L_d', 'cos_d', 'S_d', 'meff_d', 'S_r', 'eps')}
for n in np.unique(P['n']):
    f = np.load(os.path.join(eb.W_SAME['RANS'], f'fields_{n}.npz'))
    b = np.load(f'/mnt/Tank/DNS/analysis/budgets_c50b/b006/budget_{n}.npz')
    sel = P['n'] == n
    for q in ('rho_d', 'k_d', 'L_d', 'cos_d', 'S_d', 'meff_d', 'S_r'):
        d50[q][sel] = f[q].ravel()[P['idx'][sel]]
    d50['eps'][sel] = b['eps'].ravel()[P['idx'][sel]]
print('50M read', flush=True)
np.savez(out, reg=P['reg'], w=P['w'], regnames=np.array(names), blk=P['n'],
         **{'c193_' + q: v for q, v in ref.items()}, **{'rans_' + q: v for q, v in R.items()}, **d50)
print('wrote', out)
