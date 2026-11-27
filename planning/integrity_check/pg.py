"""Search extracted literature text by regex and print page number and context.

Usage: python3 -I pg.py "<file stem>" "<regex>" [context chars]
"""
import re
import sys
from pathlib import Path

here = Path(__file__).resolve().parent / "lit_txt"
stem, pat = sys.argv[1], sys.argv[2]
ctx = int(sys.argv[3]) if len(sys.argv) > 3 else 300
for f in sorted(here.glob(f"{stem}*.txt")):
    pages = f.read_text(errors="ignore").split("\f")
    for i, p in enumerate(pages, 1):
        flat = re.sub(r"\s+", " ", p)
        for m in re.finditer(pat, flat, flags=re.I):
            a, b = max(0, m.start() - ctx), min(len(flat), m.end() + ctx)
            print(f"== {f.stem} p.{i}: ...{flat[a:b]}...\n")
