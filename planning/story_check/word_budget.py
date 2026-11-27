"""Word count of manuscript/main_claude.tex per section against a budget for a 15-page paper.

Draws word_budget.pdf (bar chart, current against budget) and writes word_budget.txt (the table used in
story_check.tex). Run from any folder: python3 -I word_budget.py
Assumption: words scale with pages at fixed figure count; the budget total is current words x 15/18.
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
TEX = HERE.parent.parent / 'manuscript' / 'main_claude.tex'
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']

# share of the body word budget per block (judgement: results carry the paper; 4.1 is the longest)
SHARE = {'Abstract': None, 'Introduction': 0.15, 'Aims': 0.02, 'Methods': 0.16, 'Results lead-in': 0.03,
         '4.1 Coolant flow': 0.27, '4.2 Breakdown': 0.22, '4.3 A posteriori': 0.11, 'Conclusions': 0.04}


def words(txt):
    txt = re.sub(r'(?<!\\)%.*', '', txt)                                   # comments
    txt = re.sub(r'\\iffalse.*?\\fi', '', txt, flags=re.S)
    txt = re.sub(r'\\begin\{(equation|multline|table\*?|figure\*?|tabular)\}.*?\\end\{\1\}', ' ', txt, flags=re.S)
    txt = re.sub(r'\$[^$]*\$', ' X ', txt)                                  # inline maths = one word
    txt = re.sub(r'\\(cite|ref|label)\{[^}]*\}', ' ', txt)
    txt = re.sub(r'\\[A-Za-z]+\*?', ' ', txt)
    txt = re.sub(r'[{}\[\]~]', ' ', txt)
    return len([w for w in txt.split() if re.search(r'[A-Za-z0-9]', w)])


def main():
    src = TEX.read_text()
    abstract = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', src, re.S).group(1)
    body = src[src.index('\\section{Introduction}'):src.index('\\section*{Acknowledgments}')]
    cuts = [('Introduction', '\\section{Introduction}'), ('Aims', '\\section{Research Aims}'),
            ('Methods', '\\section{Methods}'), ('Results lead-in', '\\section{Results}'),
            ('4.1 Coolant flow', '\\subsection{Coolant Flow}'), ('4.2 Breakdown', '\\subsection{Breakdown'),
            ('4.3 A posteriori', '\\subsection{A Posteriori'), ('Conclusions', '\\section{Conclusions}')]
    pos = [body.index(c[1]) for c in cuts] + [len(body)]
    cur = {'Abstract': words(abstract)}
    for i, (name, _) in enumerate(cuts):
        cur[name] = words(body[pos[i]:pos[i + 1]])
    total_body = sum(v for k, v in cur.items() if k != 'Abstract')
    target = round(total_body * 15 / 18, -1)
    bud = {k: (200 if k == 'Abstract' else round(SHARE[k] * target, -1)) for k in cur}
    lines = [f'# words: current, budget (body budget = {target:.0f} = current body {total_body} x 15/18)']
    for k in cur:
        lines.append(f'{k:18s} {cur[k]:6d} {bud[k]:6.0f} {cur[k] - bud[k]:+6.0f}')
    (HERE / 'word_budget.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))

    names = list(cur)
    fig, ax = plt.subplots(figsize=(6.3, 2.6))
    x = range(len(names))
    ax.bar([i - 0.2 for i in x], [cur[n] for n in names], 0.4, color=OKABE_ITO[5], label='manuscript now')
    ax.bar([i + 0.2 for i in x], [bud[n] for n in names], 0.4, color=OKABE_ITO[1], hatch='//', edgecolor='k',
           linewidth=0.4, label='budget, 15 pages')
    ax.set_xticks(list(x))
    short = {'Results lead-in': 'Results\nlead-in', '4.1 Coolant flow': '4.1', '4.2 Breakdown': '4.2',
             '4.3 A posteriori': '4.3'}
    ax.set_xticklabels([short.get(n, n) for n in names], fontsize=7)
    ax.set_ylabel('words', fontsize=8)
    ax.tick_params(axis='y', labelsize=7)
    ax.legend(fontsize=7, frameon=False)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(HERE / 'word_budget.pdf')


if __name__ == '__main__':
    main()
