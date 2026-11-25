"""Shared helpers of fig_meff_split.py (F15) and fig_aeff_split.py (F16): reader of the regions_*se.txt tables, the
stacked-bar panel with total marker, full +-2 SE error bar and non-clipped y limits, and the significance report.

Panel geometry: bar centre x, bar half-width 0.21, the total marker and its error bar sit just right of the bar
(x + 0.30), the total value is printed right of the error bar at the height of the marker (boxed), segment values are
printed at the segment centre and dropped when the segment is shorter than the text. Nothing is placed outside the
axes; the y limits contain every bar, marker and error bar plus room for the panel letter.
"""
import glob
import re

import numpy as np
from matplotlib.patches import Patch

A = '/home/sk999/Documents/4dns_filmcooled/analysis/'
BW = 0.42          # bar width (data units, bar spacing is 1)
XM = 0.30          # offset of marker + error bar from the bar centre
FS = 5.0           # font size of the numbers


def pick(interim, w1_glob):
    """W1 file (reference 193M b010) if it exists, else the interim file."""
    hit = sorted(glob.glob(A + w1_glob))
    return hit[0] if hit else interim


def read_se(path):
    """{region: [total, SE, v1, SE, v2, SE, ..., fV]} from a regions*_se.txt."""
    out = {}
    for ln in open(path, encoding='utf-8'):
        if ln.startswith('#') or ln.startswith('region') or not ln.strip():
            continue
        m = re.match(r'^(.*?)\s{2,}([+-].*)$', ln.rstrip())
        out[m.group(1).strip()] = [float(v) for v in m.group(2).split()]
    return out


def fmt(v):
    return f'{v:+.2f}'.replace('-', '−')


def report(title, rungs, data, regions, terms):
    """Print value, SE, |v|/SE of every rung/region/term and list the non-significant ones (|v| < 2 SE).
    terms: [(name, index of the value in the row)]."""
    print(title)
    print('rung | region | term | value | SE | |v|/SE ; NS = not significant (|v| < 2 SE)')
    ns = []
    for lab, d in zip(rungs, data):
        for key in regions:
            r = d[key]
            for name, i in terms:
                v, se = r[i], r[i + 1]
                flag = 'NS' if abs(v) < 2 * se else ''
                print(f'{lab:9s} {key:24s} {name:9s} {v:+.3f} {se:.3f} {abs(v) / se:6.1f} {flag}')
                if flag:
                    ns.append((lab, key, name, v, se))
    print(f'\nnot significant: {len(ns)} cases')
    for c in ns:
        print('NS  %-9s %-24s %-9s %+.3f  SE %.3f' % c)
    return ns


def panel(ax, entries, colours, ylabel, letter, head=0.14, foot=0.06):
    """entries: [(x, total, SE_total, [term values])]; colours: [(colour, hatch)] per term.
    Draws stacked bars, the total marker (filled if |total| >= 2 SE, hollow otherwise), the full error bar and the
    numbers; sets y limits from the extents of everything drawn."""
    lo_all = hi_all = 0.0
    segs_all = []
    for x, tot, se, vals in entries:
        pos = neg = 0.0
        segs = []
        for v, (col, hat) in zip(vals, colours):
            bot = pos if v >= 0 else neg
            ax.bar(x, v, BW, bottom=bot, color=col, hatch=hat, edgecolor='k', linewidth=0.4)
            segs.append((v, bot + v / 2, abs(v)))
            if v >= 0:
                pos += v
            else:
                neg += v
        xm = x + XM
        ax.errorbar(xm, tot, yerr=2 * se, fmt='none', ecolor='k', capsize=1.8, lw=0.8, zorder=5)
        sig = abs(tot) >= 2 * se
        ax.plot(xm, tot, 'D', ms=3.3, mfc='k' if sig else 'white', mec='k', mew=0.8, zorder=6)
        lo_all = min(lo_all, neg, tot - 2 * se)
        hi_all = max(hi_all, pos, tot + 2 * se)
        segs_all.append((x, tot, segs))
    span = hi_all - lo_all
    ax.set_ylim(lo_all - foot * span, hi_all + head * span)
    ax.axhline(0, color='k', lw=0.6)
    fig = ax.figure
    fig.canvas.draw()
    h_pt = ax.get_window_extent().height * 72 / fig.dpi
    unit = (ax.get_ylim()[1] - ax.get_ylim()[0]) / h_pt          # data units per point
    for x, tot, segs in segs_all:
        for v, yc, h in segs:
            if h > 7.0 * unit:                                    # text height about 5 pt: drop shorter segments
                ax.text(x, yc, fmt(v), ha='center', va='center', fontsize=FS, zorder=7,
                        bbox=dict(boxstyle='square,pad=0.08', fc='white', ec='none', alpha=0.7))
        ax.text(x + XM + 0.07, tot, fmt(tot), ha='left', va='center', fontsize=FS, zorder=7,
                bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='k', lw=0.4))
    ax.set_ylabel(ylabel)
    ax.text(0.008, 0.985, letter, transform=ax.transAxes, ha='left', va='top', fontsize=7.5)


def term_handles(terms):
    """terms: [(label, colour, hatch)]."""
    return [Patch(facecolor=c, edgecolor='k', hatch=h, linewidth=0.4, label=lab) for lab, c, h in terms]
