"""Ladder coverage grid for the verification report (planning/verification/verification.tex).

Each row is a place in main_claude.tex where the RANS 50M is compared with the 193M 3DNS. Two columns give
whether the 20M LES is also compared there: in the figure or table, and in the text. Status codes:
Y compared, P partial, N not compared (new gap), K known pending, - not applicable.
Run: python3 ladder_grid.py  (writes ladder_grid.pdf beside this script).
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params  # noqa: E402

OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
STATUS = {  # colour, hatch, label
    'Y': ('#009E73', '', 'compared'),
    'P': ('#E69F00', '//', 'partial'),
    'N': ('#D55E00', 'xx', 'not compared (new gap)'),
    'K': ('#56B4E9', '..', 'known pending'),
    '-': ('#FFFFFF', '', 'not applicable'),
}

# (label, line in main_claude.tex, figure/table status, text status)
ROWS = [
    ('Abstract: effectiveness', 'L93', '-', 'Y'),
    ('Abstract: vertical flux ordering', 'L94', '-', 'N'),
    ('Abstract: eddy viscosity', 'L96', '-', 'Y'),
    ('Fig. 2 resolution (RANS 5M)', 'L392', 'K', 'Y'),
    ('Fig. 3 3D overview', 'L436', 'K', 'K'),
    ('Circulation at x/D 4', 'L440', '-', 'K'),
    ('Hole: separation, BR', 'L446', 'P', 'P'),
    ('Fig. 4 effectiveness', 'L458', 'Y', 'Y'),
    ('Fig. 5a CRVP circulation', 'L499', 'Y', 'Y'),
    ('Fig. 5a WV circulation', 'L503', 'Y', 'N'),
    ('CRVP core height', 'L501', '-', 'N'),
    ('Fig. 5b vertical flux', 'L529', 'Y', 'Y'),
    ('Flux at z/D 0.25, 0.5, 1.0', 'L535', '-', 'N'),
    ('Fig. 5c swap test', 'L540', 'P', 'P'),
    ('Fig. 5d lateral flux', 'L551', 'Y', 'P'),
    ('Upwash and inward velocity', 'L562', '-', 'N'),
    ('Plume half-width, centroid', 'L567', '-', 'N'),
    ('Mid-pitch H and Cf', 'L582', '-', 'Y'),
    ('Intermittency', 'L588', 'Y', 'P'),
    ('Fig. 7 eddy-viscosity maps', 'L669', 'Y', 'Y'),
    ('Fig. 8a log split', 'L679', 'P', 'Y'),
    ('Fig. 8b,c factor ratios', 'L682', 'N', 'N'),
    ('Cross-planes, length scale', 'L706', '-', 'N'),
    (r'Fig. 9 and $\Pr_t$ values', 'L718', 'K', 'K'),
    ('Spanwise diffusivities', 'L746', '-', 'N'),
    ('Table 2 a posteriori', 'L780', 'P', 'Y'),
    ('Fig. 10 a posteriori', 'L809', 'Y', 'Y'),
    ('Design guidance', 'L857', '-', 'Y'),
    ('Conclusions 1, 2, 4', 'L885', '-', 'Y'),
    (r'Conclusion 3 ($\Pr_t$)', 'L899', '-', 'K'),
]


def main():
    set_mpl_params(10 / 22, grid_on=False)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    n = len(ROWS)
    fig, ax = plt.subplots(figsize=(6.3, 0.22 * n + 1.0))
    fs = plt.rcParams['font.size']
    for i, (lab, line, sf, st) in enumerate(ROWS):
        y = n - 1 - i
        for j, s in enumerate((sf, st)):
            col, hatch, _ = STATUS[s]
            ax.add_patch(Rectangle((j, y), 1, 1, facecolor=col, edgecolor='k', hatch=hatch, lw=0.5))
            ax.text(j + 0.5, y + 0.5, s, ha='center', va='center', fontsize=fs - 1,
                    bbox=dict(boxstyle='square,pad=0.08', fc='white', ec='none'))
        ax.text(-0.1, y + 0.5, f'{lab} ({line})', ha='right', va='center', fontsize=fs - 1)
    ax.text(0.5, n + 0.2, 'figure or\ntable', ha='center', va='bottom', fontsize=fs - 1)
    ax.text(1.5, n + 0.2, 'text', ha='center', va='bottom', fontsize=fs - 1)
    handles = [Rectangle((0, 0), 1, 1, facecolor=c, edgecolor='k', hatch=h, lw=0.5, label=f'{k}: {l}')
               for k, (c, h, l) in STATUS.items()]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(2.15 / 2.0, 1.0), frameon=False,
              fontsize=fs - 1, handlelength=1.4)
    ax.set_xlim(0, 2)
    ax.set_ylim(0, n + 1.2)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.subplots_adjust(left=0.40, right=0.62, top=0.99, bottom=0.01)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ladder_grid.pdf')
    fig.savefig(out, bbox_inches='tight')
    print('saved', out)


if __name__ == '__main__':
    main()
