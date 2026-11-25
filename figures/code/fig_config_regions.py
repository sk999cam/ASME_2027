"""Configuration, boundary conditions and analysis regions in the centre plane, to scale (paper Fig. F1).

Usage: python3 fig_config_regions.py [OUT.pdf]
       default OUT: /home/sk999/Documents/ASME_2027/figures/fig_config_regions.pdf
       copy of this script: ASME_2027/figures/code/fig_config_regions.py
(a) whole domain: boundary conditions and the five regions above the plate.
(b) enlarged view of the hole: metering, diffuser and exit regions, coolant inlet face.
Hole outline: straight-line centre-plane outline from write_up/methods/01_configuration.tex (to be replaced by
cl_outline.npy). Region bounds: dns/python/resolution_regions.py. Downstream regions are the near film and far film.
Style: set_mpl_params, Okabe-Ito fills, no titles (CLAUDE.md figure rules).
"""
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params

OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
K, ORA, SKY, GRN, YEL, BLU, VER, PUR = OKABE_ITO
OUT = sys.argv[1] if len(sys.argv) > 1 else '/home/sk999/Documents/ASME_2027/figures/fig_config_regions.pdf'

HOLE = [(-5.446, -2.567), (-4.946, -3.433), (-2.781, -2.183), (2.36, 0.0), (-1.0, 0.0)]
METER = [(-5.446, -2.567), (-4.946, -3.433), (-2.781, -2.183), (-3.281, -1.317)]
DIFF = [(-3.281, -1.317), (-2.781, -2.183), (1.182, -0.5), (-1.866, -0.5)]
EXIT = [(-1.866, -0.5), (1.182, -0.5), (2.36, 0.0), (2.36, 0.5), (-1.0, 0.5), (-1.0, 0.0)]
X_LE, X_TE, XB, Z_FS = -1.00, 2.36, 10.25, 4.0
BUF_IN, BUF_OUT = -12.1, 22.1


def draw_geometry(ax, labels):
    """Region fills, region boundaries, buffers, walls and the coolant inlet face, in D."""
    regions = [  # (polygon, colour, label, label position)
        ([(-15, Z_FS), (25, Z_FS), (25, 10), (-15, 10)], GRN, 'freestream', (5, 7.2)),
        ([(-15, 0), (X_LE, 0), (X_LE, Z_FS), (-15, Z_FS)], SKY, 'approach layer', (-6.5, 2.0)),
        ([(X_LE, 0.5), (X_TE, 0.5), (X_TE, Z_FS), (X_LE, Z_FS)], VER, 'jet', (0.68, 2.3)),
        ([(X_TE, 0), (XB, 0), (XB, Z_FS), (X_TE, Z_FS)], BLU, 'near film', (6.3, 2.0)),
        ([(XB, 0), (25, 0), (25, Z_FS), (XB, Z_FS)], PUR, 'far film', (17.5, 2.0)),
        (METER, YEL, None, None),
        (DIFF, ORA, None, None),
        (EXIT, K, None, None),
    ]
    for poly, c, lab, pos in regions:
        ax.add_patch(Polygon(poly, closed=True, fc=c, ec='none', alpha=0.18 if c == K else 0.30))
        if lab and labels:
            ax.text(*pos, lab, ha='center', va='center', fontsize=8)
    for x in (X_LE, X_TE, XB):
        ax.plot([x, x], [0.5 if x != XB else 0, Z_FS], color='0.45', lw=0.6)
    ax.plot([-15, 25], [Z_FS, Z_FS], color='0.45', lw=0.6)
    ax.plot([-1.0, 2.36], [0.5, 0.5], color='0.45', lw=0.6)
    ax.plot([-1.866, 1.182], [-0.5, -0.5], color='0.45', lw=0.6)
    ax.plot([-1.0, -1.0], [0, 0.5], color='0.45', lw=0.6)
    ax.plot([2.36, 2.36], [0, 0.5], color='0.45', lw=0.6)
    ax.plot([-3.281, -2.781], [-1.317, -2.183], color='0.3', lw=0.8, ls='--')
    # buffers (3DNS), excluded from all data
    for x0, x1 in ((-15, BUF_IN), (BUF_OUT, 25)):
        ax.add_patch(Rectangle((x0, 0), x1 - x0, 10, fc='white', ec='none'))
        ax.add_patch(Rectangle((x0, 0), x1 - x0, 10, fc='none', ec='0.55', hatch='////', lw=0))
        if labels:
            ax.text((x0 + x1) / 2, 10.35, 'buffer', ha='center', va='bottom', fontsize=7)
    # walls, domain, coolant inlet face
    ax.plot([-15, -1.0], [0, 0], color=K, lw=1.6)
    ax.plot([2.36, 25], [0, 0], color=K, lw=1.6)
    hx, hz = zip(*HOLE)
    ax.plot([hx[0], hx[4]], [hz[0], hz[4]], color=K, lw=1.6)
    ax.plot(hx[1:4], hz[1:4], color=K, lw=1.6)
    ax.plot([hx[0], hx[1]], [hz[0], hz[1]], color=GRN, lw=2.2)
    ax.plot([-15, -15, 25, 25], [0, 10, 10, 0], color=K, lw=0.8)


def main():
    set_mpl_params(0.5, grid_on=False)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.2, 5.6), gridspec_kw=dict(height_ratios=[1.0, 0.75]))
    fs = 7
    # (a) whole domain with boundary conditions
    draw_geometry(ax, labels=True)
    ax.text(-15.6, 5, 'inlet:\nBlasius,\nRe$_\\theta$ 135-140,\nFST Tu 5.5 %', ha='right', va='center',
            fontsize=fs)
    ax.text(25.6, 7.0, 'outlet:\n$p=p_\\infty$', ha='left', va='center', fontsize=fs)
    ax.text(5, 10.35, 'top: $p=p_\\infty$ (characteristic in 3DNS)', ha='center', va='bottom', fontsize=fs)
    ax.text(17.5, -0.4, 'adiabatic no-slip plate', ha='center', va='top', fontsize=fs)
    ax.text(5, 8.8, 'span $6D$, periodic in $y$', ha='center', fontsize=fs)
    for x, s in ((X_LE, '$-1.00$'), (X_TE, '2.36'), (XB, '10.25')):
        ax.text(x, Z_FS + 0.25, s, ha='center', va='bottom', fontsize=6.5, color='0.3')
    ax.text(25.6, Z_FS, '$z/D=4$', ha='left', va='center', fontsize=6.5, color='0.3')
    ax.annotate('', (-15, -4.6), (25, -4.6), arrowprops=dict(arrowstyle='<->', lw=0.7))
    ax.text(5, -4.6, '$40D$', ha='center', va='center', fontsize=fs, bbox=dict(fc='white', ec='none', pad=1))
    ax.add_patch(Rectangle((-6.5, -3.9), 10.5, 5.5, fc='none', ec='0.3', lw=0.7, ls=':'))
    ax.text(-6.6, -3.9, '(b)', ha='right', va='bottom', fontsize=fs)
    ax.text(-23.5, 10.8, '(a)', fontsize=8)
    ax.set_xlim(-24, 31)
    ax.set_ylim(-5.4, 11.2)
    ax.set_aspect('equal')
    ax.axis('off')
    # (b) enlarged view of the hole
    draw_geometry(bx, labels=False)
    kw = dict(fontsize=7.5, va='center', arrowprops=dict(arrowstyle='-', color='0.3', lw=0.6))
    bx.annotate('metering', (-4.3, -2.6), (-4.0, -3.6), **kw)
    bx.annotate('diffuser', (-1.6, -1.1), (0.2, -1.9), **kw)
    bx.annotate('exit', (0.9, 0.2), (2.9, -0.9), **kw)
    bx.text(0.68, 1.2, 'jet', ha='center', fontsize=7.5)
    bx.text(-4.0, 1.0, 'approach layer', ha='center', fontsize=7.5)
    bx.text(3.2, 1.0, 'near film', ha='left', fontsize=7.5)
    bx.annotate('coolant inlet:\n$p_{0,c}$, $T_{0,c}$, Tu$_c$', (-5.3, -2.95), (-6.4, -1.4), color=GRN,
                fontsize=fs, va='center', arrowprops=dict(arrowstyle='-', color=GRN, lw=0.6))
    bx.annotate('metering/diffuser split,\n$3.5D$ along the axis\nfrom the breakout', (-3.0, -1.8),
                (-1.6, -3.4), fontsize=6.5, color='0.3', va='center',
                arrowprops=dict(arrowstyle='-', color='0.3', lw=0.6))
    bx.text(X_LE, 0.75, '$x/D=-1.00$', ha='right', fontsize=6.5, color='0.3')
    bx.text(X_TE, 0.75, '2.36', ha='left', fontsize=6.5, color='0.3')
    bx.text(-6.4, 1.45, '(b)', fontsize=8, va='top')
    bx.set_xlim(-6.5, 4.0)
    bx.set_ylim(-3.9, 1.6)
    bx.set_aspect('equal')
    bx.axis('off')
    fig.savefig(OUT, bbox_inches='tight')
    print('wrote', OUT)


if __name__ == '__main__':
    main()
