"""Overview of the mean flow, candidate layouts with the 20M LES added (previews only, nothing is written to figures/).

Usage: python3 fig_overview_les_options.py A     (option A strip; the existing two-view figure figures/fig_overview_3d.png is stacked above it)
       python3 fig_overview_les_options.py B     (option B: three 3D views, 193M, RANS 50M, 20M LES, each with its three insets)
Output: planning/review/57_overview_optionA.png, 57_overview_optionA_strip.png, 57_overview_optionB.png

Uses fig_overview_3d.py unchanged (scene(), data(), find_cores(), trace_lines(), classify_cores(), nest_save()). The 20M LES enters through the cache
files of the LES (read only): overview_3d/les20_q_omx_grid.npz (Q, omega_x on the n4 grid), mean_features/les20_mean_grid.npz (U, V, W, omega_x on the
n4 grid), and the wall eta of the plate-vertex cache of fig_eta_x ('20M LES', gridded as fig_eta_wall.py). The LES has no fine cross-planes at
x/D 2, 4, 10: its planes are slices of the n4 grid (0.05 D) interpolated linearly onto the (y, z) grid of the 193M planes (0.02 D). The LES
is treated like the RANS (no smoothing of the plotted fields; the 193M is smoothed as in fig_overview_3d.py).
"""
import io
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
import matplotlib
matplotlib.use('Agg')
import fig_overview_3d as O  # noqa: E402

REV = '/home/sk999/Documents/ASME_2027/planning/review/'
AS = '/mnt/Elements/DNS/analysis/asme2027_figs/'
SCR = '/tmp/claude-1001/-home-sk999-Documents-ASME-2027/f4c92ceb-6df8-4254-96c3-6552d15d0ecf/scratchpad/'
EXISTING = '/home/sk999/Documents/ASME_2027/figures/fig_overview_3d.png'
D, UINF = O.D, O.UINF
LAB = {1: '193M 3DNS', 0: 'RANS 50M', 2: '20M LES'}


def les_planes():
    """Slices of the LES mean grid at x/D 2, 4, 10 (omega_x D/U_inf, V, W in m/s) on the (y, z) grid of the 193M fine planes."""
    from scipy.interpolate import RegularGridInterpolator as RGI
    P = np.load(AS + 'mean_features/les20_planes_b.npz')
    if all(('om_x%d' % xp) in P.files for xp in (2, 4, 10)):          # fine planes (0.02 D) of the LES, same format as rans50_planes.npz
        return {xp: (P['om_x%d' % xp], P['V_x%d' % xp], P['W_x%d' % xp], P['y'], P['z']) for xp in (2, 4, 10)}
    M = np.load(AS + 'mean_features/les20_mean_grid.npz')           # fallback: slices of the 0.05 D grid
    n3 = np.load(O.N3)
    y, z = n3['y'], n3['z']
    Y, Z = np.meshgrid(y, z, indexing='ij')
    pl = {}
    for xp in (2, 4, 10):
        i = int(np.argmin(np.abs(M['x'] - xp)))
        assert abs(M['x'][i] - xp) < 1e-6
        f = lambda a: RGI((M['y'], M['z']), np.nan_to_num(np.asarray(a[i], float)), bounds_error=False, fill_value=None)(np.c_[Y.ravel(), Z.ravel()]).reshape(Y.shape)  # noqa: E731
        pl[xp] = (f(M['omx_n']), f(M['V']), f(M['W']), y, z)
    return pl


def inject_les(d):
    """Replace the RANS entries of the shared data dict by the LES ones (scene(0, cfg) then draws the LES). Returns a restore function."""
    import fig_eta_x as F
    import fig_eta_wall as FW
    saved = dict(R=d['R'], eta=d['eta'], planes={xp: d['planes'][xp]['R'] for xp in d['planes']})
    LG = np.load(AS + 'overview_3d/les20_q_omx_grid.npz')
    ex, ey, ee = F.plate_points('20M LES')
    eta_les = np.ma.filled(FW.grid_eta(np.asarray(ex), np.asarray(ey), np.asarray(ee)), np.nan)
    d['R'] = LG
    d['eta'] = (eta_les, d['eta'][1])
    pl = les_planes()
    for xp in (2, 4, 10):
        d['planes'][xp]['R'] = pl[xp]
    i20 = int(np.argmin(np.abs(d['V']['x'] - 20.0)))
    d['planes'][20]['R'] = (np.array(LG['omx_n'][i20]), None, None, d['V']['y'], d['V']['z'])
    d['planes_les'] = pl

    def restore():
        d['R'], d['eta'] = saved['R'], saved['eta']
        for xp, v in saved['planes'].items():
            d['planes'][xp]['R'] = v
    return restore


def inject_les_as193(d):
    """Same data, put into the slots of the 193M so that scene(1, cfg) draws the LES with the smoothing the 193M gets (sm_q for the Q
    isosurface, sm for the plotted omega_x and (V, W) of the planes, sm_iso for the surface colour). Returns a restore function."""
    import fig_eta_x as F
    import fig_eta_wall as FW
    saved = dict(V=d['V'], eta=d['eta'], planes={xp: d['planes'][xp]['D'] for xp in d['planes']})
    LG = np.load(AS + 'overview_3d/les20_q_omx_grid.npz')
    ex, ey, ee = F.plate_points('20M LES')
    eta_les = np.ma.filled(FW.grid_eta(np.asarray(ex), np.asarray(ey), np.asarray(ee)), np.nan)
    V = d['V']
    d['V'] = dict(x=V['x'], y=V['y'], z=V['z'], Q_n=np.array(LG['Q_n']), omx=np.array(LG['omx_n']) * UINF / D)
    d['eta'] = (d['eta'][0], eta_les)
    pl = les_planes()
    for xp in (2, 4, 10):
        d['planes'][xp]['D'] = pl[xp]
    i20 = int(np.argmin(np.abs(V['x'] - 20.0)))
    d['planes'][20]['D'] = (np.array(LG['omx_n'][i20]), None, None, V['y'], V['z'])

    def restore():
        d['V'], d['eta'] = saved['V'], saved['eta']
        for xp, v in saved['planes'].items():
            d['planes'][xp]['D'] = v
    return restore


def config(wpr):
    from functions import set_mpl_params
    c = dict(O.CFG)
    K = c['ss_two']
    c.update(ss=K, win=(2400 * K, 1300 * K), dpi=c['dpi_two'], final_w=int(round(600 * wpr)), two_bars=True, lab_below=0.6, tick_sp=0.3, clear=0.08,
             cap_clear=0.15)
    set_mpl_params(O.SF_TWO * 7.5 / wpr, grid_on=False)
    return c


def render(which, c, tag):
    """scene() of one rung (0 RANS, 1 193M, 2 LES), cached on disk (the render takes minutes)."""
    f = SCR + 'ovl_scene_%s.pkl' % tag
    if os.path.exists(f) and not os.environ.get('OVL_RERENDER'):
        return pickle.load(open(f, 'rb'))
    d = O.data()
    if which == 2:                                    # LES with the smoothing of the 193M (default) or raw like the RANS (OVL_LES_RAW=1)
        raw = bool(os.environ.get('OVL_LES_RAW'))
        restore = (inject_les if raw else inject_les_as193)(d)
        try:
            r = O.scene(0 if raw else 1, c)
        finally:
            restore()
    else:
        r = O.scene(which, c)
    pickle.dump(r, open(f, 'wb'))
    return r


TAG = {0: 'rans', 1: '193', 2: 'les'}


def plane_src(which, xp):
    d = O.data()
    if which == 2:
        return O.les_planes_cache()[xp] if hasattr(O, 'les_planes_cache') else _LP()[xp]
    return d['planes'][xp]['R' if which == 0 else 'D']


_LPC = {}


def _LP():
    if not _LPC:
        _LPC.update(les_planes())
    return _LPC


def inset_data(which, xp, c, zmax):
    """Same as fig_overview_3d._inset_plane for any rung: omega_x (193M and LES smoothed with sm / sm_uv), streamline paths, cores."""
    from scipy.ndimage import gaussian_filter
    o, v, w, y1, z1 = plane_src(which, xp)
    smooth = which != 0
    zi = z1 <= zmax + 1e-9
    z1, o, v, w = z1[zi], np.array(o[:, zi], float), np.array(v[:, zi], float), np.array(w[:, zi], float)
    cs = O.find_cores(o, y1, z1, c['core_thr'], c['core_sm'] if smooth else 0, zmax)
    osm = gaussian_filter(np.nan_to_num(o), c['core_sm']) if smooth else np.nan_to_num(o)
    for side in (-1, 1):
        if any(q[0] * side > 0 and abs(q[0]) < 1.3 and q[2] * (-side) > 0 for q in cs):
            continue
        iy, iz = np.where((y1[:, None] * side > 0) & (np.abs(y1[:, None]) < 1.3) & (z1[None, :] >= 0.15))
        if len(iy):
            k_ = np.argmax((-side) * osm[iy, iz])
            if (-side) * osm[iy[k_], iz[k_]] >= c['core_thr']:
                cs.append((y1[iy[k_]], z1[iz[k_]], osm[iy[k_], iz[k_]]))
    if smooth:
        m = np.isnan(o)
        o = gaussian_filter(np.nan_to_num(o), c['sm'])
        o[m] = np.nan
        v, w = gaussian_filter(np.nan_to_num(v), c['sm_uv']), gaussian_filter(np.nan_to_num(w), c['sm_uv'])
    seeds = []
    for yc, zc, _ in cs:
        for r in (0.08, 0.16):
            for t in np.linspace(0, 2 * np.pi, 6, endpoint=False):
                seeds.append((yc + r * np.cos(t), zc + 0.5 * r * np.sin(t)))
    seeds += [(a_, b_) for a_ in np.arange(-2.9, 3.0, 0.4) for b_ in np.arange(0.08, zmax, 0.16)]
    paths = O.trace_lines(y1, z1, v, w, seeds, c['inset_cell'] * (c['cell2'] if xp == 2 else 1.0), zmax, c['smin'], ds=0.02, lmax=2.5, lmin=0.2)
    return o, paths, cs, y1, z1


def wv_items(cs, c):
    cs = [q for q in cs if q[1] <= 1.0]
    return [it for it in O.classify_cores(cs) if it[0] == 'wall' and abs(it[1][2]) >= c['label_thr']]


def draw_inset(ia, which, xp, c, zmax, label=True, fs=10):
    from matplotlib.colors import Normalize
    o_, paths_, cs_, y_, z_ = inset_data(which, xp, c, zmax)
    ia.pcolormesh(y_, z_, np.nan_to_num(o_).T, cmap='bwr', norm=Normalize(-c['omax'], c['omax']), shading='auto', rasterized=True)
    for pth in paths_:
        ia.plot(pth[:, 0], pth[:, 1], color='#202020', lw=0.6)
        m_ = len(pth) // 2
        dv = pth[min(m_ + 1, len(pth) - 1)] - pth[max(m_ - 1, 0)]
        dv = dv / (np.hypot(*dv) + 1e-12)
        ia.annotate('', xy=pth[m_] + 0.025 * dv, xytext=pth[m_] - 0.025 * dv, arrowprops=dict(arrowstyle='-|>', lw=0, color='#202020', mutation_scale=5))
    ia.set_xlim(3, -3)                                  # y/D increases to the left: looking downstream, as in the 3D view
    ia.set_ylim(0, zmax)
    ia.set_aspect('equal')
    out = []
    if label:
        for kind, q, side in wv_items(cs_, c):
            out.append((q[0], q[1]))
            ia.annotate('WV', xy=(q[0], q[1]), xytext=(q[0], zmax * 0.97), fontsize=fs, ha='center', va='top', annotation_clip=False,
                        bbox=dict(fc='w', ec='none', pad=0.5, alpha=0.75), arrowprops=dict(arrowstyle='-|>', lw=0.6, color='k', mutation_scale=5))
    return out


def strip(c, zmax=1.5):
    """Option A strip: rows 193M, RANS 50M, 20M LES; columns x/D 2, 4, 10. Axis labels on the outer panels only."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib import cm
    W = 7.5
    wc, gx, gy = 1.90, 0.26, 0.12
    xl0 = 1.12                                         # left margin: row name, z/D label, tick labels
    hc = wc * zmax / 6.0
    top, bot = 0.26, 0.52
    H = top + 3 * hc + 2 * gy + bot
    fig = plt.figure(figsize=(W, H))
    rows = [(1, '193M 3DNS'), (0, 'RANS 50M'), (2, '20M LES')]
    for ir, (wh, nm) in enumerate(rows):
        for ic, xp in enumerate((2, 4, 10)):
            x0 = xl0 + ic * (wc + gx)
            y0 = bot + (2 - ir) * (hc + gy)
            ax = fig.add_axes([x0 / W, y0 / H, wc / W, hc / H])
            draw_inset(ax, wh, xp, c, zmax, fs=10)
            ax.set_xticks([-3, 0, 3])
            ax.set_yticks([0, 1])
            ax.tick_params(length=2, pad=1.5)
            if ir < 2:
                ax.set_xticklabels([])
            else:
                for lab_t, ha_t in zip(ax.get_xticklabels(), ('right', 'center', 'left')):
                    lab_t.set_ha(ha_t)
            if ic > 0:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel(r'$z/D$', labelpad=1)
            if ir == 2 and ic == 1:
                ax.set_xlabel(r'$y/D$', labelpad=1)
            if ir == 0:
                ax.text(0.5, 1.04, r'$x/D=%d$' % xp, transform=ax.transAxes, ha='center', va='bottom')
            if ic == 0:
                fig.text(0.03 / W, (y0 + 0.5 * hc) / H, nm, ha='left', va='center')
    print('strip size %.2f x %.2f in; each inset %.2f x %.2f in, z/D 0-%.1f' % (W, H, wc, hc, zmax))
    return fig, H


def option_a(c):
    import matplotlib.pyplot as plt
    from PIL import Image
    fig, H = strip(c)
    b = io.BytesIO()
    fig.savefig(b, dpi=300, format='png', facecolor='white')
    plt.close(fig)
    st = Image.open(b).convert('RGB')
    st.save(REV + '57_overview_optionA_strip.png', dpi=(300, 300))
    ex = Image.open(EXISTING).convert('RGB')                      # existing two-view figure, unchanged (600 dpi, 6.5 in wide)
    exs = ex.resize((int(ex.size[0] / 2), int(ex.size[1] / 2)))   # 300 dpi
    canvas = Image.new('RGB', (st.size[0], exs.size[1] + st.size[1] + 30), 'white')
    canvas.paste(exs, ((st.size[0] - exs.size[0]) // 2, 0))
    canvas.paste(st, (0, exs.size[1] + 30))
    canvas = canvas.resize((canvas.size[0] * 2 // 3, canvas.size[1] * 2 // 3), Image.LANCZOS)    # 200 dpi preview
    canvas.save(REV + '57_overview_optionA.png', dpi=(200, 200))
    print('option A preview', canvas.size, 'at 200 dpi; strip %.2f in high; existing figure %.2f x %.2f in' % (H, ex.size[0] / 600, ex.size[1] / 600))


def compose3(res, c, order=(1, 0, 2)):
    """fig_overview_3d.compose for three panels (copy of its layout logic, no change of the algorithm): res[which] = (image, anchors, cores) for
    which = 0 RANS, 1 193M, 2 LES; order = top-to-bottom. Every caption is placed at the top left of its panel."""
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import Normalize
    W = 5.6
    K = c.get('ss', 1)
    rl = [(res[w][0], res[w][1] / K, res[w][2]) for w in order]
    cols = [np.where((im < 250).any((0, 2)))[0] / K for im, _, _ in rl]
    c0, c1 = min(q.min() for q in cols) - 10, max(q.max() for q in cols) + 10
    crops = []
    for im, an, _ in rl:
        rr = np.where((im < 250).any((1, 2)))[0] / K
        crops.append((rr.min() - c['pad_top'], max(rr.max() + c['pad_bot'], an[0][1] + c['pad_hs'])))
    hpx = max(b - a for a, b in crops)
    wi_ = c['img_w']
    hp = wi_ * hpx / (c1 - c0)
    wi_in = c['inset_w']
    hi_in = wi_in * c['zmax_in'] / 6.0
    lay = []
    for k_, (im_, an_, _) in enumerate(rl):
        r1o = crops[k_][0] + hpx
        A = [((a[0] - c0) / (c1 - c0) * wi_, (r1o - a[1]) / hpx * hp) for a in (an_[8], an_[9], an_[10])]
        lb, la = np.polyfit([a_[0] for a_ in A], [a_[1] for a_ in A], 1)
        row = []
        for i_ in range(3):
            xl = c['xs_in'][i_]
            top = -c['clr'] if i_ == 0 else la + lb * xl - c['clr']
            row.append((xl, top - hi_in))
        lay.append(row)
    hr = max(0.05, c['tick_sp'] - min(p_[1] for l_ in lay for p_ in l_))
    mt = 0.15
    Hp = hp + hr + mt
    figs = [plt.figure(figsize=(W, Hp)) for _ in rl]
    fs = plt.rcParams['font.size']
    caps = {1: '(a) 193M 3DNS', 0: '(b) RANS 50M', 2: '(c) 20M LES'}
    for k, (fig, (im, an, _cores)) in enumerate(zip(figs, rl)):
        which = order[k]
        ax = fig.add_axes([0.0, hr / Hp, wi_ / W, hp / Hp])
        pd_ = 400
        im = np.pad(im, ((pd_ * K, pd_ * K), (0, 0), (0, 0)), constant_values=255)
        an = an + np.array([0, pd_])
        r0 = crops[k][0] + pd_
        r1 = r0 + hpx
        ax.imshow(im[int(round(r0 * K)):int(round(r1 * K)), int(round(c0 * K)):int(round(c1 * K))], extent=(c0, c1, r1, r0))
        ax.set_xlim(c0, c1)
        ax.set_ylim(r1, r0)
        ax.axis('off')
        ax.text(c['tag_xy'][0], c['tag_xy'][1], caps[which], transform=ax.transAxes, va='top', ha=c['tag_ha'])
        hs, c1p, c2p, t2, t4, t10, t20, tcl, s2p, s1p, s10p = an
        ax.annotate('HV', xy=hs, xytext=(hs[0] - 150, hs[1] + 120), arrowprops=dict(arrowstyle='->', lw=1.0, color='k'), va='center', ha='center')
        H = Hp
        pbot = hr
        for i_in, (xp, zmax, anchor) in enumerate(((2, c['zmax_in'], s2p), (4, c['zmax_in'], s1p), (10, c['zmax_in'], s10p))):
            wi = c['inset_w']
            hi = wi * zmax / 6.0
            x0in, yb_ = lay[k][i_in]
            ibot = pbot + yb_
            ia = fig.add_axes([x0in / W, ibot / H, wi / W, hi / H])
            o_, paths_, cs_, y_, z_ = inset_data(which, xp, c, zmax)
            ia.pcolormesh(y_, z_, np.nan_to_num(o_).T, cmap='bwr', norm=Normalize(-c['omax'], c['omax']), shading='auto', rasterized=True)
            for pth in paths_:
                ia.plot(pth[:, 0], pth[:, 1], color='#202020', lw=0.6)
                m_ = len(pth) // 2
                dv = pth[min(m_ + 1, len(pth) - 1)] - pth[max(m_ - 1, 0)]
                dv = dv / (np.hypot(*dv) + 1e-12)
                ia.annotate('', xy=pth[m_] + 0.025 * dv, xytext=pth[m_] - 0.025 * dv,
                            arrowprops=dict(arrowstyle='-|>', lw=0, color='#202020', mutation_scale=5))
            ia.set_xlim(3, -3)
            ia.set_ylim(0, zmax)
            ia.set_aspect('equal')
            ia.set_xticks([-3, 0, 3])
            ia.set_yticks([0, 2])
            ia.tick_params(labelsize=fs, length=2, pad=1.5)
            for lab_t, ha_t in zip(ia.get_xticklabels(), ('right', 'center', 'left')):
                lab_t.set_ha(ha_t)
            if i_in == 0:
                ia.set_xlabel(r'$y/D$', fontsize=fs, labelpad=1)
                ia.set_ylabel(r'$z/D$', fontsize=fs, labelpad=1)
            xa = 0.0
            ax_in = (anchor[0] - c0) / (c1 - c0) * wi_
            cx = min(max(ax_in, x0in + 0.05), x0in + wi - 0.05)
            ax.annotate('', xy=anchor, xycoords='data', xytext=(cx / W, (ibot + hi) / H), textcoords='figure fraction',
                        arrowprops=dict(arrowstyle='-', lw=0.6, color='0.25', shrinkA=0, shrinkB=0))
            kk = dict(fontsize=fs, ha='center', va='top', arrowprops=dict(arrowstyle='-|>', lw=0.6, color='k', mutation_scale=6))
            items = wv_items(cs_, c)
            xt = np.array([it[1][0] for it in items] + [xa], float)
            wd = np.array([0.9 for it in items] + [0.3])
            for _ in range(300):
                order_ = np.argsort(xt)
                for a_, b_ in zip(order_[:-1], order_[1:]):
                    dd = xt[b_] - xt[a_]
                    fa, fb = a_ == len(xt) - 1, b_ == len(xt) - 1
                    gp = 0.5 * (wd[a_] + wd[b_]) + c['lab_pad']
                    if dd < gp:
                        sh = gp - dd + 1e-3
                        if fa:
                            xt[b_] += sh
                        elif fb:
                            xt[a_] -= sh
                        else:
                            xt[a_] -= 0.5 * sh
                            xt[b_] += 0.5 * sh
                xt[:-1] = np.clip(xt[:-1], -c['lab_clip'], c['lab_clip'])
            tgt = np.array([it[1][0] for it in items], float)
            for _ in range(50):
                o2 = np.argsort(tgt)
                for a_, b_ in zip(o2[:-1], o2[1:]):
                    if tgt[b_] - tgt[a_] < 0.7:
                        sh = 0.7 - (tgt[b_] - tgt[a_]) + 1e-3
                        tgt[a_] -= 0.5 * sh
                        tgt[b_] += 0.5 * sh
            for i, (kind, q, side_) in enumerate(items):
                ia.annotate('WV', xy=(tgt[i], q[1]), xytext=(xt[i], -c['lab_below']), textcoords='data', xycoords='data', annotation_clip=False, **kk)
        mid = 0.5 * (c1p + c2p)
        lx, ly = mid[0] + c['cvp_dx'], mid[1] + c['cvp_dy']
        ax.text(lx, ly, 'CRVP', ha='center', va='bottom')
        for cp_ in (c1p, c2p):
            ax.annotate('', xy=cp_, xytext=(lx, ly + 12), arrowprops=dict(arrowstyle='->', lw=1.0, color='k', shrinkA=0, shrinkB=1))
        ax.text(t2[0] - 14, t2[1] + 6, r'$x/D=2$', ha='right', va='top')
        for t, s in ((t4, r'$x/D=4$'), (t20, r'$x/D=20$')):
            ax.text(t[0], t[1] - 15, s, ha='center', va='bottom')
        ax.text(t10[0] - 12, t10[1] - 15, r'$x/D=10$', ha='right', va='bottom')
        ax.text(tcl[0] - 10, tcl[1], r'$y=0$', ha='right', va='center')
    strips = []
    Hs = 0.62
    for x_, nrm_, cmp_, lab_ in ((0.71, Normalize(-c['omax'], c['omax']), 'bwr', r'$\bar{\omega}_x D/U_\infty$'), (0.44, Normalize(0, O.VMAX), 'jet', r'$\eta$')):
        fs_ = plt.figure(figsize=(W, Hs))
        cx_ = fs_.add_axes([x_, 0.40 / Hs, 0.22, 0.12 / Hs])
        fs_.colorbar(cm.ScalarMappable(norm=nrm_, cmap=cmp_), cax=cx_, orientation='horizontal', label=lab_)
        strips.append(fs_)
    O.nest_save(figs + strips, c, caption=None)


def option_b(wpr):
    c = config(wpr)
    c['out'] = REV + '57_overview_optionB.png'
    res = {w: render(w, c, TAG[w]) for w in (0, 1, 2)}
    compose3(res, c)


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'A'
    c = config(6.5)
    if mode == 'A':
        option_a(c)
    elif mode == 'B':
        option_b(float(sys.argv[2]) if len(sys.argv) > 2 else 4.6)
    elif mode == 'scenes':
        for w in (1, 0, 2):
            im, an, cs = render(w, c, TAG[w])
            print(LAB[w], im.shape)
