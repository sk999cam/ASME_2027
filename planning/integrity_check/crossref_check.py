"""Compare cited bib entries with their Crossref records (DOI lookups only; no other data sent).

Usage: python3 -I crossref_check.py > crossref_check.txt
"""
import json
import re
import time
import urllib.request
from pathlib import Path

MS = Path(__file__).resolve().parents[2] / "manuscript"
tex = (MS / "main_claude.tex").read_text()
cited = set()
for line in tex.splitlines():
    if line.lstrip().startswith("%"):
        continue
    for m in re.finditer(r"\\cite\{([^}]*)\}", line):
        cited.update(k.strip() for k in m.group(1).split(","))

bib = (MS / "references.bib").read_text() + "\n" + (MS / "references_claude.bib").read_text()
entries = {}
for chunk in re.split(r"\n(?=@)", bib):
    m = re.match(r"\s*@(\w+)\{([^,]+),", chunk)
    if not m:
        continue
    fields = dict(
        (k.lower(), v.strip())
        for k, v in re.findall(r"\n\s*(\w+)\s*=\s*\{(.*?)\},?\s*(?=\n)", chunk, flags=re.S)
    )
    entries[m.group(2)] = (m.group(1), fields)


def norm(s):
    return re.sub(r"[^a-z0-9]", "", re.sub(r"\\[a-z]+|[{}]", "", s.lower()))


for key in sorted(cited):
    typ, f = entries[key]
    doi = f.get("doi")
    if not doi:
        print(f"{key}\tNO_DOI\t{typ}\t{f.get('title','')}")
        continue
    url = "https://api.crossref.org/works/" + urllib.request.quote(doi)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "bibcheck/1.0"})
        msg = json.load(urllib.request.urlopen(req, timeout=20))["message"]
    except Exception as e:  # noqa: BLE001
        print(f"{key}\tLOOKUP_FAILED\t{doi}\t{e}")
        continue
    ct = (msg.get("title") or [""])[0]
    cy = msg.get("issued", {}).get("date-parts", [[None]])[0][0]
    pub = msg.get("published-print", msg.get("published-online", {})).get("date-parts", [[None]])[0][0]
    cv, cp = msg.get("volume"), msg.get("page") or msg.get("article-number")
    cauth = [a.get("family", "") for a in msg.get("author", [])]
    issues = []
    if norm(ct) != norm(f.get("title", "")):
        issues.append(f"title: bib '{f.get('title','')}' vs crossref '{ct}'")
    if str(cy) != f.get("year") and str(pub) != f.get("year"):
        issues.append(f"year: bib {f.get('year')} vs crossref issued {cy}/print {pub}")
    if f.get("volume") and cv and f.get("volume") != cv:
        issues.append(f"volume: bib {f.get('volume')} vs crossref {cv}")
    if f.get("pages") and cp and norm(f["pages"]) != norm(cp):
        issues.append(f"pages: bib {f.get('pages')} vs crossref {cp}")
    if not f.get("pages") and cp:
        issues.append(f"pages missing in bib; crossref {cp}")
    if not f.get("volume") and cv and typ == "article":
        issues.append(f"volume missing in bib; crossref {cv}")
    bfirst = norm(f.get("author", "").split(",")[0])
    if cauth and bfirst != norm(cauth[0]):
        issues.append(f"first author: bib {f.get('author','')[:30]} vs crossref {cauth[0]}")
    print(f"{key}\t{'OK' if not issues else 'DIFF'}\t{doi}\t" + " | ".join(issues))
    time.sleep(0.3)
