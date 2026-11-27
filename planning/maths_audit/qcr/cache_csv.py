#!/usr/bin/env python3
"""Read one rans50 grad CSV ONCE (streaming, pyarrow) and cache every column as float32 .npy (grad-file row order).
usage: python3 -B cache_csv.py TAG     (TAG = base | qcr | qcr_prt04)
Input read-only: /mnt/Elements/DNS/analysis/rans50_exports/rans50_TAG_cells_grad.csv
Output: /mnt/Elements/DNS/analysis/asme2027_figs/qcr/TAG/<col>.npy  (+ DONE marker)
"""
import os
import sys
import numpy as np
import pyarrow as pa
import pyarrow.csv as pc

NAMES = ['u', 'v', 'w', 'dUdx', 'dUdy', 'dUdz', 'dVdx', 'dVdy', 'dVdz', 'dWdx', 'dWdy', 'dWdz', 'Smag', 'T',
         'dTdx', 'dTdy', 'dTdz', 'rho', 'mu', 'k', 'omega', 'mut', 'wd', 'gamma', 'V', 'X', 'Y', 'Z']
N = 53120567
tag = sys.argv[1]
src = f'/mnt/Elements/DNS/analysis/rans50_exports/rans50_{tag}_cells_grad.csv'
out = f'/mnt/Elements/DNS/analysis/asme2027_figs/qcr/{tag}'
os.makedirs(out, exist_ok=True)
mm = {n: np.lib.format.open_memmap(f'{out}/{n}.npy', mode='w+', dtype=np.float32, shape=(N,)) for n in NAMES}
ro = pc.ReadOptions(block_size=1 << 26, column_names=NAMES, skip_rows=1)
co = pc.ConvertOptions(column_types={n: pa.float64() for n in NAMES})
rd = pc.open_csv(src, read_options=ro, convert_options=co)
i = 0
for b in rd:
    m = b.num_rows
    for n in NAMES:
        mm[n][i:i + m] = b.column(n).to_numpy().astype(np.float32)
    i += m
    print(i, flush=True)
assert i == N, (i, N)
for n in NAMES:
    mm[n].flush()
open(f'{out}/DONE', 'w').write(f'{i}\n')
print('done', i)
