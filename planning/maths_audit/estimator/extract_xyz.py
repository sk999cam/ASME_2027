#!/usr/bin/env python3
"""Coordinates of the point set of pts.npz (50M grid only; no 193M read). Writes xyz/D to OUT.npz and checks the order against pts.npz.
usage: PYTHONDONTWRITEBYTECODE=1 python3 -B extract_xyz.py pts.npz OUT.npz"""
import os
import sys
import numpy as np
os.environ['MEFF_REGIONS'] = 'pp'
os.environ['MEFF_SFX_C193b'] = '_b010'
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/sk999/Documents/4dns_filmcooled/analysis/meff_c50b_c100b')
import errorbars as eb  # noqa: E402
names, P = eb.points(eb.CASES['C50b'], 2)
d = np.load(sys.argv[1])
assert np.array_equal(P['reg'], d['reg']) and np.array_equal(P['n'], d['blk']) and np.allclose(P['w'], d['w']), 'order mismatch'
np.savez(sys.argv[2], xyz=P['xyz'] / eb.D if hasattr(eb, 'D') else P['xyz'] / 0.00775, idx=P['idx'])
print('ok', len(P['n']))
