#!/usr/bin/env python3
"""The 13 classes of the mesh (meshing/coarsning topology): one bold line per class, coloured by its
direction (i streamwise, j spanwise, k normal to the plate / along the hole), over a faint wireframe (every 6th grid line).

    python3 plot_mesh_classes.py [MESH.cgns [OUT.pdf]]
    default MESH: /mnt/Elements/DNS/meshes/rans_20M_Lm25b.cgns
    default OUT:  /home/sk999/Documents/ASME_2027/figures/fig_mesh_classes.pdf (paper Fig. F2)
    copy of this script: ASME_2027/figures/code/fig_mesh_classes.py

(a) plate plane, all 23 blocks: the plate-plane and hole cross-section classes.
(b) centre plane through the hole: the wall-normal class and the two classes along the hole.
Each class's line is picked in PICK (zone, line, target point); block boundaries in black.
Style: set_mpl_params, Okabe-Ito colours with a line style per direction, no titles (CLAUDE.md figure rules).
"""
import os, sys

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params                 # noqa: E402

OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
OUT_DEFAULT = '/home/sk999/Documents/ASME_2027/figures/fig_mesh_classes.pdf'
MESH_DEFAULT = '/mnt/Elements/DNS/meshes/rans_20M_Lm25b.cgns'

ROOT = '/home/sk999/Documents/4dns_filmcooled'          # absolute, so the copy in figures/code runs as well
sys.path.insert(0, os.path.join(ROOT, 'meshing', 'coarsning', 'python'))
from coarsen import config as C                      # noqa: E402
from coarsen.cgnsio import read_topology             # noqa: E402
from coarsen.topology import classes                 # noqa: E402

D = 7.75e-3
STRIDE = 6
GCOL = {0: OKABE_ITO[5], 1: OKABE_ITO[3], 2: OKABE_ITO[6]}   # i blue, j green, k vermillion
GLS = {0: '-', 1: '--', 2: '-.'}                               # and a line style per direction
SIDE = ('col-dom-17', 'col-dom-1', 'col-dom-3', 'col-dom-9', 'col-dom-12', 'col-dom-15',
        'hole-core', 'hole-LE', 'hole-TE')


def main():
    mesh = sys.argv[1] if len(sys.argv) > 1 else MESH_DEFAULT
    out = sys.argv[2] if len(sys.argv) > 2 else OUT_DEFAULT
    set_mpl_params(0.5, grid_on=False)
    plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
    f = h5py.File(mesh, 'r')
    topo = read_topology(mesh)
    src = read_topology(os.path.join(ROOT, 'meshing', 'stage9_export.cgns'))
    cl, bp = classes(src)
    X = {zn: np.stack([f[f'Base/{zn}/GridCoordinates/Coordinate{q}/ data'][...] for q in 'XYZ'], -1)
         for zn in topo}                            # (nk, nj, ni, 3)
    kpl = X['col-dom-1'].shape[0] - 1               # plate plane (0-based)
    # coarse breakpoints = the coarse file's patch ends; the same segments, in the same order
    cbp = {}
    for zn, z in topo.items():
        cbp[zn] = []
        for d in range(3):
            s = {1, z.dims[d]}
            for p in z.bcs + z.conns:
                s.update(p.span(d))
            cbp[zn].append(sorted(s))
    where = {}
    for c in cl:
        for m in c.members:
            k = bp[m.seg.zone][m.seg.d].index(m.seg.a)
            where[(m.seg.zone, m.seg.d, k)] = c
    assert all(len(cbp[zn][d]) == len(bp[zn][d]) for zn in src for d in range(3)), 'breakpoints differ'

    fig = plt.figure(figsize=(11, 6.4))
    ax1 = fig.add_axes([0.05, 0.58, 0.93, 0.38]); ax2 = fig.add_axes([0.05, 0.08, 0.62, 0.42])
    axl = fig.add_axes([0.70, 0.08, 0.29, 0.42]); axl.axis('off'); axl.set_xlim(0, 1); axl.set_ylim(0, 1)

    def strided(P, axis, u, v):
        """Every STRIDE-th index along `axis`, counted from the end with the smaller coordinate in the
        direction the index runs (blocks such as col-dom-18 run backwards), so faint lines line up
        across block faces; both end lines always included."""
        n = P.shape[axis]
        a, b = np.take(P, 0, axis), np.take(P, n - 1, axis)
        dvec = (b - a).reshape(-1, P.shape[-1]).mean(0)[[u, v]]
        rev = dvec[np.argmax(np.abs(dvec))] < 0
        qs = set(range(0, n, STRIDE)) | {0, n - 1}
        return sorted(n - 1 - q for q in qs) if rev else sorted(qs)

    def faint(ax, P, u, v):
        for q in strided(P, 0, u, v):
            ax.plot(P[q, :, u] / D, P[q, :, v] / D, c='0.78', lw=0.3)
        for q in strided(P, 1, u, v):
            ax.plot(P[:, q, u] / D, P[:, q, v] / D, c='0.78', lw=0.3)

    best = {}                                        # class id -> (length, panel, line)

    # where each class's bold line is drawn (class number 1-13): (zone, line index across, target x/D, y/D or z/D)
    # zone None = any member; line None = every line, the one passing nearest the target is taken
    PICK = {1: ('hole-TE', 0, (2.4, 0.0)),      # hole TE lip (wall line of the TE O-block)
            9: ('hole-jp', 0, (0.7, 0.6)),      # hole left-hand lip (wall line of the +y O-block)
            11: ('hole-LE', None, (-1.0, 0.0)), # O-block wall-normal at the LE
            2: ('col-dom-4', None, (0.6, -1.6)), # collar radial, bottom (-y) side
            8: (None, None, (-10.0, 0.0)), 10: (None, None, (-3.3, 0.0)),   # streamwise on y = 0
            5: (None, None, (6.5, 0.0)), 7: (None, None, (17.0, 0.0)),
            4: ('col-dom-7', None, (-3.3, -1.7)), 6: ('col-dom-8', None, (-3.3, 1.7)),  # spanwise i lines, both sides
            # (not upstream of x/D -4.6: there col-dom-18 has i and j swapped)
            3: (None, None, (2.7, 1.0)),
            12: ('hole-core', 'mid', (0.0, -1.0)), 13: ('hole-core', 'mid', (-4.0, -2.4))}  # hole axis

    def candidates(panel, zn, P, dI, dJ, u, v):
        """P (nJ, nI, 3); lines along I (direction dI) and along J (direction dJ)."""
        for dA, axis in ((dI, 1), (dJ, 0)):
            ends = cbp[zn][dA]
            n_other = P.shape[0] if axis == 1 else P.shape[1]
            for s, (a0, a1) in enumerate(zip(ends[:-1], ends[1:])):
                c = where[(zn, dA, s)]
                if (C.class_type(c) in ('fs_k', 'hole_ax')) != (panel == 2):
                    continue                         # k classes on (b), in-plane classes on (a)
                zone, line, tgt = PICK[c.id + 1]
                if zone is not None and zn != zone:
                    continue
                qs = range(n_other) if line is None else [n_other // 2 if line == 'mid' else line]
                for q in qs:
                    L = P[q, a0 - 1:a1] if axis == 1 else P[a0 - 1:a1, q]
                    L = L[:, [u, v]] / D
                    score = -np.hypot(L[:, 0] - tgt[0], L[:, 1] - tgt[1]).min()
                    if score > best.get(c.id, (-np.inf,))[0]:
                        best[c.id] = (score, panel, L, dA)   # dA: index direction in its block

    def edges(ax, P, u, v):
        for L in (P[0], P[-1], P[:, 0], P[:, -1]):
            ax.plot(L[:, u] / D, L[:, v] / D, c='k', lw=0.9, zorder=3)

    for zn in topo:
        P = X[zn][kpl]                               # (nj, ni, 3)
        faint(ax1, P, 0, 1)
        edges(ax1, P, 0, 1)
        candidates(1, zn, P, 0, 1, 0, 1)
    for zn in SIDE:
        Q = X[zn][:, :, X[zn].shape[2] // 2]         # (nk, nj, 3)
        faint(ax2, Q, 0, 2)
        edges(ax2, Q, 0, 2)
        candidates(2, zn, Q, 1, 2, 0, 2)
    axes = {1: ax1, 2: ax2}
    placed = {1: [], 2: []}                          # label positions, kept >= 0.6 D apart
    for c in cl:
        _, panel, L, g = best[c.id]
        # display convention (not the blocks' own index labels): i streamwise, j spanwise,
        # k normal to the plate or along the hole; from the drawn line's extent
        if panel == 2:
            g = 2
        else:
            g = 0 if np.ptp(L[:, 0]) >= np.ptp(L[:, 1]) else 1
        ax = axes[panel]
        ax.plot(L[:, 0], L[:, 1], c=GCOL[g], ls=GLS[g], lw=2.6, solid_capstyle='butt',
                dash_capstyle='butt', zorder=5)
        lim = {1: ((-15, 25), (-3, 3)), 2: ((-8, 6), (-3.5, 2.2))}[panel]
        V = L[(L[:, 0] > lim[0][0]) & (L[:, 0] < lim[0][1]) & (L[:, 1] > lim[1][0]) & (L[:, 1] < lim[1][1] - 0.3)]
        cum = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(V, axis=0), axis=1))])
        for t in (0.5, 0.3, 0.7, 0.2, 0.8, 0.1, 0.9):
            m = V[np.searchsorted(cum, t * cum[-1]).clip(0, len(V) - 1)]
            if all(np.hypot(*(m - o)) >= 0.6 for o in placed[panel]):
                break
        placed[panel].append(m)
        box = dict(boxstyle='circle,pad=0.15', fc='white', ec=GCOL[g], lw=1)
        if cum[-1] < 0.8:                            # short line (O-block): label to the side, with a leader
            ax.annotate(str(c.id + 1), xy=m, xytext=(m[0] - 1.6, m[1] - 1.2), fontsize=7, ha='center',
                        va='center', zorder=6, bbox=box,
                        arrowprops=dict(arrowstyle='-', color=GCOL[g], lw=0.9))
        else:
            ax.text(m[0], m[1], str(c.id + 1), fontsize=7, ha='center', va='center', zorder=6, bbox=box)
    ax1.set_aspect('equal'); ax1.set_xlim(-15, 25); ax1.set_ylim(-3, 3)
    ax1.set_xlabel('x/D'); ax1.set_ylabel('y/D')
    ax1.text(0.0, 1.02, '(a)', transform=ax1.transAxes, ha='left', va='bottom')
    ax2.set_aspect('equal'); ax2.set_xlim(-8, 6); ax2.set_ylim(-3.5, 2.2)
    ax2.set_xlabel('x/D'); ax2.set_ylabel('z/D')
    ax2.text(0.0, 1.02, '(b)', transform=ax2.transAxes, ha='left', va='bottom')
    y = 0.85
    for g, txt in ((0, 'i: streamwise'), (1, 'j: spanwise'), (2, 'k: normal to the plate,\n    along the hole')):
        axl.plot([0, 0.12], [y] * 2, c=GCOL[g], ls=GLS[g], lw=2.6, dash_capstyle='butt')
        axl.text(0.16, y, txt, va='center')
        y -= 0.13
    axl.text(0, y - 0.03, 'grey: every 6th grid line', va='top', color='0.4')
    fig.savefig(out)
    print('wrote', out, {c.id + 1: round(-best[c.id][0], 2) for c in cl})


if __name__ == '__main__':
    main()
