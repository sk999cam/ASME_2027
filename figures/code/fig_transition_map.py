"""Transition mechanism map on the HPT suction surface against Re_x (paper Fig. F1, Introduction).

Reused, not redrawn (user decision 2026-11-24): Fig. 1 of Kar et al. (2026, GPPS-TC-2026-0005), itself reproduced
from Clark et al. (2024). The figure is cut out of page 2 of the published PDF as a vector crop, so it keeps the
original styling; it is the one paper figure not drawn with set_mpl_params. Cite both papers in the caption.

Usage: python3 fig_transition_map.py
       writes /home/sk999/Documents/ASME_2027/figures/fig_transition_map.pdf
Needs ghostscript (gs).
"""
import subprocess

SRC = '/home/sk999/Documents/4dns_filmcooled/write_up/literature/papers/reference/Kar et al. 2026.pdf'
OUT = '/home/sk999/Documents/ASME_2027/figures/fig_transition_map.pdf'
PAGE = 2
# crop box in PDF points (origin bottom left), A4 page 595.28 x 841.89: the figure without its caption
X0, Y0, X1, Y1 = 178.0, 422.0, 428.0, 607.0

subprocess.run(['gs', '-q', '-o', OUT, '-sDEVICE=pdfwrite', f'-dFirstPage={PAGE}', f'-dLastPage={PAGE}',
                '-dUseCropBox',
                '-c', f'[/CropBox [{X0} {Y0} {X1} {Y1}] /PAGES pdfmark',
                '-f', SRC], check=True)
print('wrote', OUT)
