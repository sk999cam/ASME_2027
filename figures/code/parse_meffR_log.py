"""Parse estimator_modeC.log (planning/maths_audit/estimator) into two machine-readable tables for fig_meff_split.py.

Usage: python3 parse_meffR_log.py [log_path]    (default log: estimator_modeC.log; far film x/D 10-20: estimator_modeC_x20.log)
Writes (figures/code/data/):
  meffR_modeC.txt : per region the 4-player Shapley terms (section B), the 5-player L_eff and align terms and their
                    order ranges [min, max] over the 120 orders.
  meffR_se.txt    : SE of the total per region (section F, delete-one spatial-block jackknife). To replace the SE
                    (time-batch jackknife), overwrite this file with the same columns and change the '# kind:' line.
Regions: jet x/D 0-3, film 3-10, film 10-25 (plate regions only).
"""
import os
import re
import sys

# optional argument: path of the estimator log (default: estimator_modeC.log); the far-film region keeps the log label 'film 10-25'
LOG = sys.argv[1] if len(sys.argv) > 1 else '/home/sk999/Documents/ASME_2027/planning/maths_audit/estimator/estimator_modeC.log'
LOGNAME = os.path.basename(LOG)
DATA = '/home/sk999/Documents/ASME_2027/figures/code/data/'
REG = ['jet 0-3', 'film 3-10', 'film 10-25']
NUM = r'([+-]\d+\.\d+)'

lines = open(LOG, encoding='utf-8').read().splitlines()


def block(start_pat, stop_pat=None):
    i = next(n for n, l in enumerate(lines) if re.search(start_pat, l))
    out = []
    for l in lines[i + 1:]:
        if stop_pat and re.search(stop_pat, l):
            break
        out.append(l)
    return out


def row(blk, name):
    return next(l for l in blk if re.match(r'^' + re.escape(name) + r'\s', l))


# section B: 5-player (rho, energy, L_eff, align, S-wt, L_al, sum-tot) and 4-player (rho, energy, L_al, S-wt)
B = block(r'^B\. Exact additive split', r'^5 players: min/max')
# order ranges, 5 players
R5 = block(r'^5 players: min/max', r'^4 players \(L_al')
# section F: total +- SE
F = block(r'^F\. Headline 4-player', None)

res = {}
for r in REG:
    l = row(B, r)
    left, mid, right = l.split('|')
    tot = float(left.split()[-1])
    f5 = [float(v) for v in re.findall(NUM, mid)]            # rho energy L_eff align S-wt L_al
    f4 = [float(v) for v in re.findall(NUM, right)]          # rho energy L_al S-wt
    rl = row(R5, r)
    rng = [(float(a), float(b)) for a, b in re.findall(r'\[' + NUM + ',' + NUM + r'\]', rl)]  # rho en L al S
    fl = row(F, r)
    m = re.search(NUM + r'\+-(\d+\.\d+)', fl)
    assert abs(float(m.group(1)) - tot) < 1e-9
    res[r] = dict(total=tot, se=float(m.group(2)), rho=f4[0], energy=f4[1], L_al=f4[2], S=f4[3],
                  L_eff=f5[2], align=f5[3], L_eff_rng=rng[2], align_rng=rng[3])

with open(DATA + 'meffR_modeC.txt', 'w', encoding='utf-8') as f:
    f.write(f'# parsed from {LOGNAME} (mode C: wall nodes excluded, z/D >= 0.3); ln(mu_R,RANS / mu_R,193)\n')
    f.write('# 4-player Shapley: density, energy, L_al, S-wt; secondary 5-player L_eff, align with order ranges\n')
    f.write('region total density energy L_al Swt L_eff L_eff_min L_eff_max align align_min align_max\n')
    for r in REG:
        d = res[r]
        f.write(f"{r.replace(' ', '_')} {d['total']:+.3f} {d['rho']:+.3f} {d['energy']:+.3f} {d['L_al']:+.3f} "
                f"{d['S']:+.3f} {d['L_eff']:+.3f} {d['L_eff_rng'][0]:+.3f} {d['L_eff_rng'][1]:+.3f} "
                f"{d['align']:+.3f} {d['align_rng'][0]:+.3f} {d['align_rng'][1]:+.3f}\n")
with open(DATA + 'meffR_se.txt', 'w', encoding='utf-8') as f:
    f.write('# kind: spatial-block jackknife\n')
    f.write(f'# SE of the total ln(mu_R,RANS / mu_R,193); source {LOGNAME} section F\n')
    f.write('region SE\n')
    for r in REG:
        f.write(f"{r.replace(' ', '_')} {res[r]['se']:.3f}\n")
for r in REG:
    print(r, res[r])
