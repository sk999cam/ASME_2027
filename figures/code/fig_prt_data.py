"""F19: turbulent Prandtl number and directional diffusivity ratios against the data (paper Fig. F19, INTERIM).

Usage: python3 fig_prt_data.py
Output: /home/sk999/Documents/ASME_2027/figures/fig_prt_data.pdf

(a) Like-for-like Pr_t,z = nu_z/alpha_z (the eps_m/eps_Ty form of Kohli & Bogard 2005, exp of the volume-averaged
    ln) on the centre plane (|y|/D <= 0.06), x/D 1-10, in the z/D bands 0-0.25, 0.25-0.5, 0.5-1, 1-2, for the
    50M 3DNS and the 20M LES; error bars +-2 jackknife SE (6 batches), horizontal bars the band. The RANS value 0.9
    is the line. LES points are shifted by +0.04 in z/D to separate them; the values are at the band centres.
    Kohli & Bogard (2005) data are NOT drawn: no digitised file of them exists in the repo (only the PDF), and
    nothing is invented. Their Pr_t of about 2 is quoted in the text only.
(b) ln(nu_y/nu_z) (momentum, filled) and ln(alpha_y/alpha_z) (heat, open) per region, +-2 SE, 50M 3DNS and
    20M LES. Grey band: Kaszeta & Simon (2000) momentum ratio 1.5-4.0, i.e. ln 0.41-1.39 (their y and z swapped
    against the project axes). Not drawn: metering (0.10 / 0.05 of the volume passes the guards) and approach layer
    (0.02 / none), listed as not a result in analysis/directional/README.md.

INTERIM: no 193M value (the directional analysis exists for the 50M and the 20M LES only, stat3 b006). The regions
are the point-by-point regions at stride 2 plus the centre-plane points; the guards (0.05 U/D, 0.05 dT0/D) are those of
the README. Data: /home/sk999/Documents/4dns_filmcooled/analysis/directional/results/dir.txt (produced by
dir_table.py, 50M 3DNS and 20M LES at stat3 b006).
Style: set_mpl_params, Okabe-Ito rung colours (50M #D55E00, 20M LES #E69F00, RANS #009E73), no title.
"""
import re
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler

sys.path.insert(0, '/home/sk999/Documents/3DNS_FP')
from functions import set_mpl_params

SRC = '/home/sk999/Documents/4dns_filmcooled/analysis/directional/results/dir.txt'
OUT = '/home/sk999/Documents/ASME_2027/figures/fig_prt_data.pdf'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
C50, CLES, CRANS = '#D55E00', '#E69F00', '#009E73'


def parse():
    """{rung: {'reg': {(region, quantity): (value, SE)}, 'cp': [(z0, z1, Pr, SE)]}}"""
    out, cur = {}, None
    for ln in open(SRC, encoding='utf-8'):
        if ln.startswith('## '):
            cur = out.setdefault(ln[3:].strip(), {'reg': {}, 'cp': []})
            continue
        m = re.match(r'^cp z/D ([\d.]+)-([\d.]+)\s+Pr_t,z\s+([\d.]+) \(SE ([\d.]+)\)', ln)
        if m and cur is not None:
            cur['cp'].append(tuple(float(v) for v in m.groups()))
            continue
        m = re.match(r'^(.{11}) (ln\(\S+\))\s+([+-][\d.]+)\s+([\d.]+)', ln)
        if m and cur is not None:
            cur['reg'][(m.group(1).strip(), m.group(2))] = (float(m.group(3)), float(m.group(4)))
    return out


set_mpl_params(0.4, grid_on=True)
plt.rcParams['axes.prop_cycle'] = cycler(color=OKABE_ITO)
D = parse()
fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw={'width_ratios': [1, 1.25]})

# ---- (a)
for rung, col, ls, mk, off, lab in (('50M 3DNS', C50, ':', 's', 0.0, '50M 3DNS'),
                                    ('20M LES', CLES, '-.', 'o', 0.04, '20M LES')):
    cp = D[rung]['cp']
    zc = np.array([(c[0] + c[1]) / 2 for c in cp]) + off
    hw = np.array([(c[1] - c[0]) / 2 for c in cp])
    pr = np.array([c[2] for c in cp])
    se = np.array([c[3] for c in cp])
    a.errorbar(zc, pr, yerr=2 * se, xerr=hw, fmt=mk, color=col, ls=ls, capsize=2, lw=1.0, label=lab)
a.axhline(0.9, color=CRANS, ls='-', label=r'RANS $\mathrm{Pr}_t=0.9$')
a.set_xlim(0, 2.05)
a.set_ylim(0, 1.1)
a.set_xlabel(r'$z/D$')
a.set_ylabel(r'$\mathrm{Pr}_{t,z}=\nu_z/\alpha_z$')
a.text(0.03, 0.97, '(a)', transform=a.transAxes, ha='left', va='top')
a.legend(loc='lower right', fontsize=8)

# ---- (b)
REG = [('diffuser', 'diffuser'), ('exit', 'exit'), ('jet 0-3', 'jet'), ('film 3-10', 'near\nfilm'),
       ('film 10-25', 'far\nfilm')]
b.axhspan(np.log(1.5), np.log(4.0), color='0.8', zorder=0, label=r'Kaszeta & Simon (2000), $\nu$: 1.5–4.0')
b.axhline(0, color='k', lw=0.6)
slots = [('50M 3DNS', 'ln(nu_y/nu_z)', C50, 's', -0.27, True, r'50M, $\nu_y/\nu_z$'),
         ('50M 3DNS', 'ln(al_y/al_z)', C50, 's', -0.09, False, r'50M, $\alpha_y/\alpha_z$'),
         ('20M LES', 'ln(nu_y/nu_z)', CLES, 'o', 0.09, True, r'20M LES, $\nu_y/\nu_z$'),
         ('20M LES', 'ln(al_y/al_z)', CLES, 'o', 0.27, False, r'20M LES, $\alpha_y/\alpha_z$')]
for rung, q, col, mk, off, filled, lab in slots:
    xs, vs, es = [], [], []
    for i, (key, _) in enumerate(REG):
        v, s = D[rung]['reg'][(key, q)]
        xs.append(i + off); vs.append(v); es.append(2 * s)
    b.errorbar(xs, vs, yerr=es, fmt=mk, color=col, mfc=col if filled else 'white', ls='', capsize=2, lw=1.0, label=lab)
b.set_xticks(range(len(REG)))
b.set_xticklabels([r[1] for r in REG])
b.set_xlim(-0.55, len(REG) - 0.45)
b.set_ylabel(r'$\ln(\nu_y/\nu_z)$, $\ln(\alpha_y/\alpha_z)$')
b.text(0.97, 0.97, '(b)', transform=b.transAxes, ha='right', va='top')
h, l = b.get_legend_handles_labels()
b.legend(h[1:] + h[:1], l[1:] + l[:1], loc='upper left', fontsize=6.5, ncol=1, bbox_to_anchor=(0.0, 0.93))
b.set_ylim(-1.0, 2.2)
fig.tight_layout()
fig.savefig(OUT, bbox_inches='tight')
print('wrote', OUT)
