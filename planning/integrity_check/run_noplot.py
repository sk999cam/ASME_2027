"""Run a figure script for its printed output only: every savefig is replaced by a no-op, so no file is written.

Usage: python3 -B run_noplot.py <script.py> [args]
Only for scripts whose sole file output is savefig (checked by grep before use).
"""
import runpy
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.figure  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


def _noop(*a, **k):
    print('[savefig suppressed]', a[1] if len(a) > 1 else (a[0] if a and isinstance(a[0], str) else ''))


matplotlib.figure.Figure.savefig = _noop
plt.savefig = _noop
script = sys.argv[1]
sys.argv = sys.argv[1:]
sys.path.insert(0, '/home/sk999/Documents/ASME_2027/figures/code')
runpy.run_path(script, run_name='__main__')
