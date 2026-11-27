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
from matplotlib.patches import ConnectionPatch, Polygon, Rectangle
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
        ([(-15, Z_FS), (25, Z_FS), (25, 10), (-15, 10)], GRN, 'freestream', (5.0, 7.0)),
        ([(-15, 0), (X_LE, 0), (X_LE, Z_FS), (-15, Z_FS)], SKY, 'approach layer', (-6.5, 2.0)),
        ([(X_LE, 0.5), (X_TE, 0.5), (X_TE, Z_FS), (X_LE, Z_FS)], VER, 'jet', (0.68, 2.25)),
        ([(X_TE, 0), (XB, 0), (XB, Z_FS), (X_TE, Z_FS)], BLU, 'near film', (6.3, 2.0)),
        ([(XB, 0), (25, 0), (25, Z_FS), (XB, Z_FS)], PUR, 'far film', (17.5, 2.0)),
        (METER, YEL, None, None),
        (DIFF, ORA, None, None),
        (EXIT, K, None, None),
    ]
    for poly, c, lab, pos in regions:
        ax.add_patch(Polygon(poly, closed=True, fc=c, ec='none', alpha=0.18 if c == K else 0.30))
        if lab and labels:
            ax.text(*pos, lab, ha='center', va='center', fontsize=plt.rcParams['font.size'])
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
            ax.text((x0 + x1) / 2, 10.35, 'buffer', ha='center', va='bottom', fontsize=plt.rcParams['font.size'])
    # walls, domain, coolant inlet face
    ax.plot([-15, -1.0], [0, 0], color=K, lw=1.6)
    ax.plot([2.36, 25], [0, 0], color=K, lw=1.6)
    hx, hz = zip(*HOLE)
    ax.plot([hx[0], hx[4]], [hz[0], hz[4]], color=K, lw=1.6)
    ax.plot(hx[1:4], hz[1:4], color=K, lw=1.6)
    ax.plot([hx[0], hx[1]], [hz[0], hz[1]], color=GRN, lw=2.2)
    ax.plot([-15, -15, 25, 25], [0, 10, 10, 0], color=K, lw=0.8)


def main():
    # Drawn at the printed size, 6.75 in = 0.9 text width (figure* with width=0.9\textwidth, scale 1); body text 10 pt.
    # Side by side: whole domain (a) left, small tight zoom of the hole (b) right, both with true aspect.
    set_mpl_params(10 / 22, grid_on=False)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    FW = 6.75
    XL, XR_ = -19.0, 29.5
    YL, YH = -7.8, 12.4
    BX0, BX1, BZ0, BZ1 = -5.8, 4.0, -5.7, 0.6     # extent of the zoom = zoom box in (a)
    zoom_w = 1.8
    gap = 0.2
    main_w = FW - zoom_w - gap
    m_sc = main_w / (XR_ - XL)         # inches per D in the main view
    z_sc = zoom_w / (BX1 - BX0)        # inches per D in the zoom
    zoom_h = (BZ1 - BZ0) * z_sc
    fh = (YH - YL) * m_sc
    fw = FW
    fig = plt.figure(figsize=(fw, fh))
    ax = fig.add_axes([0, 0, main_w / fw, 1.0])
    bx = fig.add_axes([(main_w + gap) / fw, (fh - zoom_h) / 2 / fh, zoom_w / fw, zoom_h / fh])
    fs = plt.rcParams['font.size']
    fs_s = fs - 1
    # (a) whole domain with boundary conditions
    draw_geometry(ax, labels=True)
    # boundary-condition notes (inlet, outlet, top, plate, span) removed: they are in the paper text
    for x, s, ha in ((X_LE, '$x/D=-1.00$', 'right'), (X_TE, '2.36', 'left'), (XB, '10.25', 'center')):
        ax.text(x, Z_FS + 0.15, s, ha=ha, va='bottom', fontsize=fs, color='0.3')
    ax.text(25.6, Z_FS, '$z/D=4$', ha='left', va='center', fontsize=fs, color='0.3')
    ax.annotate('', (-15, -7.0), (25, -7.0), arrowprops=dict(arrowstyle='<->', lw=0.7))
    ax.text(5, -7.0, '$40D$', ha='center', va='center', fontsize=fs, bbox=dict(fc='white', ec='none', pad=1))
    ax.add_patch(Rectangle((BX0, BZ0), BX1 - BX0, BZ1 - BZ0, fc='none', ec='0.3', lw=0.7, ls=':'))
    ax.text(BX0 - 0.2, BZ0, '(b)', ha='right', va='bottom', fontsize=fs)
    ax.text(-18.8, 11.8, '(a)', fontsize=fs, va='top')
    ax.set_xlim(XL, XR_)
    ax.set_ylim(YL, YH)
    ax.set_aspect('equal')
    ax.axis('off')
    # (b) enlarged view of the hole
    draw_geometry(bx, labels=False)
    kw = dict(fontsize=fs, va='center', arrowprops=dict(arrowstyle='-', color='0.3', lw=0.6))
    # only what is not labelled in (a): metering, diffuser, exit, coolant inlet, split note
    bx.annotate('metering', (-4.0, -3.0), (-3.0, -4.0), ha='center', **kw)
    bx.annotate('diffuser', (-0.5, -1.4), (1.6, -2.0), ha='center', **kw)
    bx.annotate('exit', (1.0, 0.0), (3.2, -0.9), ha='center', **kw)
    bx.annotate('coolant inlet', (-5.3, -3.3), (-5.7, -5.2), color=GRN, ha='left',
                fontsize=fs, va='center', arrowprops=dict(arrowstyle='-', color=GRN, lw=0.6))
    bx.annotate('split at $3.5D$\nalong the axis', (-3.0, -1.9),
                (2.0, -3.9), fontsize=fs_s, color='0.3', va='center', ha='center',
                arrowprops=dict(arrowstyle='-', color='0.3', lw=0.6))
    bx.text(BX0 + 0.1, -0.75, '(b)', fontsize=fs, va='center')
    bx.set_xlim(BX0, BX1)
    bx.set_ylim(BZ0, BZ1)
    bx.set_aspect('equal')
    bx.axis('off')
    fig.savefig(OUT)
    print('wrote', OUT)


if __name__ == '__main__':
    main()
