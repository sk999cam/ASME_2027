"""Bar chart of the integrity-check counts (citations A3 and reported numbers B).

Usage: python3 fig_counts.py
Output: planning/integrity_check/fig_counts.pdf (beside integrity_check.tex)
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params  # noqa: E402

OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_counts.pdf')

# counts: match, rounding, mismatch, no source / source not stored, not checkable
CATS = ['match / supported', 'rounding', 'mismatch / not supported', 'no source / not stored', 'not checkable']
A3 = [48, 0, 6, 15, 0]
B = [165, 5, 6, 10, 5]
COLS = [OKABE_ITO[3], OKABE_ITO[2], OKABE_ITO[6], OKABE_ITO[1], '0.6']
HATCH = ['', '//', 'xx', '..', '--']

set_mpl_params(10 / 22, grid_on=True)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
fig, ax = plt.subplots(figsize=(5.5, 2.4))
y = np.array([1.0, 0.0])
left = np.zeros(2)
for i, c in enumerate(CATS):
    v = np.array([A3[i], B[i]])
    ax.barh(y, v, left=left, color=COLS[i], hatch=HATCH[i], edgecolor='k', linewidth=0.5, height=0.55, label=c)
    for yy, ll, vv in zip(y, left, v):
        if vv >= 20:
            ax.text(ll + vv / 2, yy, str(vv), ha='center', va='center', color='k')
    left += v
ax.set_yticks(y)
ax.set_yticklabels([f'A3 claims (n={sum(A3)})', f'B numbers (n={sum(B)})'])
ax.set_xlabel('count')
ax.set_xlim(0, 200)
ax.legend(loc='upper center', bbox_to_anchor=(0.45, -0.35), ncol=3, frameon=False)
fig.tight_layout()
fig.savefig(OUT)
print('saved', OUT)
