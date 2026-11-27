"""F8 fig_flow_features: instantaneous centre-plane temperature of the 193M 3DNS with the jet shear-layer roll-up (Results 3).

Usage: python3 fig_flow_features.py      (writes /home/sk999/Documents/ASME_2027/figures/fig_flow_features.png)

(a) instantaneous 193M centre plane y = 0, full stored extent x/D -2..10, z/D -1..3: temperature
    (T0,fs - T)/(T0,fs - T0,c) (jet) with the solid in light grey, the wall/hole/diffuser outline, the Kelvin-Helmholtz roll-up of
    the windward shear layer and the cross-plane stations x/D = 2, 4, 10 of the overview figure (fig_overview_3d).
(b) zoom of the near-hole region x/D -1..6 (frame in (a)).
Frame: record 1999, step 75670, t U/D = 249.09 (not printed on the figure).
The stored frame ends at x/D = 10 (z/D -1..3); nothing is padded or extrapolated.

Data: /mnt/Elements/DNS/analysis/flow_features_c193b/kh_frame_c193b.npz (producer resample_n3n4.py of 4dns_filmcooled).
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from cycler import cycler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params  # noqa: E402

OUT = '/home/sk999/Documents/ASME_2027/figures/fig_flow_features.png'
FRAME = '/mnt/Elements/DNS/analysis/flow_features_c193b/kh_frame_c193b_x20.npz'   # instantaneous centre-plane frame (one constant)
XCAP = 20.0
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
PLANES = (2, 4, 10)
ZLIM = (-1.0, 3.0)
ZOOM = (-1.0, 6.0)
SOLID = '0.82'


def field(ax, x, z, eta, ol, xlim, ylim):
    im = ax.pcolormesh(x, z, eta.T, cmap='jet', vmin=0, vmax=1, shading='nearest', rasterized=True)
    ax.set_facecolor(SOLID)
    ax.plot(ol[:, 0], ol[:, 1], 'k-', lw=0.8, zorder=3)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect('equal')
    ax.set_ylabel(r'$z/D$')
    return im


def main():
    K = np.load(FRAME, allow_pickle=True)
    x, z, T, ol = K['x'], K['z'], K['T'], K['outline']
    eta = (300.0 - T) / (300.0 - 199.02)
    xmax = min(float(x.max()), XCAP)
    set_mpl_params(10.0 / 22.0, grid_on=False)     # font 10 pt = body size; drawn at printed width (3.56 in), include scale 1
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    W = 3.56
    L, R = 0.50, 0.10
    wa = W - L - R
    sc = wa / (xmax - float(x.min()))                # inch per unit of x/D (equal aspect)
    ha, hb = (ZLIM[1] - ZLIM[0]) * sc, 3.0 * wa / (ZOOM[1] - ZOOM[0])
    yc, hc = 0.42, 0.12
    yb = yc + hc + 0.55
    ya = yb + hb + 0.52
    H = ya + ha + 0.05
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([L / W, ya / H, wa / W, ha / H])
    im = field(ax, x, z, eta, ol, (x.min(), xmax), ZLIM)
    ax.set_xlabel(r'$x/D$', labelpad=1)
    ax.set_yticks([0, 2])
    bz = fig.add_axes([L / W, yb / H, wa / W, hb / H])
    field(bz, x, z, eta, ol, ZOOM, (-1.0, 2.0))
    bz.set_xlabel(r'$x/D$', labelpad=1)
    ax.text(0.01, 0.95, '(a)', transform=ax.transAxes, va='top', color='w')
    bz.text(0.01, 0.95, '(b)', transform=bz.transAxes, va='top', color='w')
    # frame of (b) on (a)
    ax.add_patch(Rectangle((ZOOM[0], -1.0), ZOOM[1] - ZOOM[0], 3.0, fill=False, ec='w', lw=0.8, ls='--', zorder=4))
    # cross-plane stations
    for p in PLANES:
        xt = min(p, xmax - 0.02)
        ax.plot([xt, xt], [2.15, 3.0], 'w-', lw=1.0, zorder=4)
        ax.text(xt + 0.35, 2.95, '%d' % p, color='w', ha='left', va='top', zorder=5)
    ax.text(12.0, 2.95, 'cross planes', color='w', ha='left', va='top', zorder=5)
    # annotations on (b)
    kw = dict(arrowstyle='->', lw=0.8, color='w')
    bz.annotate('Kelvin-Helmholtz\nroll-up', xy=(0.9, 0.55), xytext=(2.6, -0.2), color='w',
                ha='left', va='center', zorder=6, arrowprops=kw)
    bz.annotate('', xy=(1.9, 0.75), xytext=(2.6, 0.2), zorder=6, arrowprops=kw)
    for p in (2, 4):
        bz.plot([p, p], [1.6, 2.0], 'w-', lw=1.0, zorder=4)
        bz.text(p, 1.55, r'$x/D=%d$' % p, color='w', ha='center', va='top', zorder=5)
    cax = fig.add_axes([L / W, yc / H, wa / W, hc / H])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    cb.set_label(r'$(T_{0,fs}-T)/(T_{0,fs}-T_{0,c})$', labelpad=1)
    fig.savefig(OUT, dpi=600)
    print('t =', float(K['t']), 's, tU/D =', float(K['tU_D']), 'record', int(K['record']), 'step', int(K['step']))


if __name__ == '__main__':
    main()
