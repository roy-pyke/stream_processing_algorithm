#!/usr/bin/env python3
"""Static checks of the LaTeX sources (no TeX installation needed), plus a page check of the built PDF.

  - every \\NumXxx macro used in paper/*.tex is defined in paper/numbers.tex and is not "??"
  - every \\ref / \\eqref target exists; references into the other document use the prefixes of the
    xr package: A-<label> (or \\apx{<label>}) in paper.tex points to appendix.tex, P-<label> (or
    \\pref{<label>}) in appendix.tex points to paper.tex
  - every \\cite key exists in paper/references.bib
  - \\begin / \\end environments and braces are balanced, $ signs are paired
  - every \\input and \\includegraphics file exists
  - double-anonymous review: no link to the authors' own repository
  - if paper/paper.pdf exists (scripts/build_paper.sh): the body ends within 12 pages (references may
    start at the top of page 13)
Exit status 1 if anything is wrong.
"""
import os
import re
import shutil
import subprocess
import sys

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper")
bad = 0


def problem(*msg):
    global bad
    print(*msg)
    bad += 1


read = lambda n: open(os.path.join(root, n)).read()
docs = {"paper": ["paper.tex", "results.tex"], "appendix": ["appendix.tex"]}
texts = {n: read(n) for files in docs.values() for n in files}
strip = lambda t: re.sub(r"(?<!\\)%.*", "", t)
body = {d: "\n".join(strip(texts[n]) for n in files) for d, files in docs.items()}

# number macros
defs = dict(re.findall(r"\\newcommand\{\\(Num[A-Za-z]+)\}\{(.*)\}", read("numbers.tex")))
used = sorted(set(re.findall(r"\\(Num[A-Za-z]+)", "\n".join(body.values()))))
for u in used:
    if u not in defs:
        problem("UNDEFINED macro:", u)
    elif "??" in defs[u] or "nan" in defs[u].lower():
        problem("macro without value:", u, "=", defs[u])

# references
labels = {d: set(re.findall(r"\\label\{([^}]+)\}", t)) for d, t in body.items()}
for d, t in body.items():
    other = "appendix" if d == "paper" else "paper"
    pre = "A-" if d == "paper" else "P-"
    refs = re.findall(r"\\(?:eq)?ref\{([^}]+)\}", t)
    refs += [pre + r for r in re.findall(r"\\(?:apx|pref)\{([^}]+)\}", t)]
    for r in sorted(set(refs)):
        if "#" in r:          # macro definitions such as \apx
            continue
        if r.startswith(pre):
            if r[len(pre):] not in labels[other]:
                problem(f"{d}: undefined cross-document reference: {r}")
        elif r not in labels[d]:
            problem(f"{d}: undefined reference: {r}")

# citations
bibkeys = set(re.findall(r"@\w+\{([^,\s]+),", read("references.bib")))
cited = set()
for c in re.findall(r"\\cite(?:\[[^\]]*\])?\{([^}]+)\}", "\n".join(body.values())):
    cited.update(x.strip() for x in c.split(","))
for c in sorted(cited - bibkeys):
    problem("undefined citation:", c)
for c in sorted(bibkeys - cited):
    print("note: bibliography entry never cited:", c)

# structure
for name, t in texts.items():
    stack = []
    for m in re.finditer(r"\\(begin|end)\{([^}]+)\}", strip(t)):
        if m.group(1) == "begin":
            stack.append(m.group(2))
        elif not stack or stack.pop() != m.group(2):
            problem(f"{name}: unbalanced environment near offset {m.start()}: {m.group(0)}")
    if stack:
        problem(f"{name}: unclosed environments: {stack}")
    b = strip(t).replace("\\{", "").replace("\\}", "")
    if b.count("{") != b.count("}"):
        problem(f"{name}: unbalanced braces: {b.count('{')} open, {b.count('}')} close")
    if len(re.findall(r"(?<!\\)\$", b)) % 2:
        problem(f"{name}: odd number of $ signs")
    if re.search(r"roy-pyke|stream_processing_algorithm", t):
        problem(f"{name}: link to the authors' repository (double-anonymous review)")
allt = "\n".join(body.values())
for f in re.findall(r"\\input\{([^}]+)\}", allt):
    if not os.path.exists(os.path.join(root, f if f.endswith(".tex") else f + ".tex")):
        problem("missing \\input file:", f)
for f in re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", allt):
    if not os.path.exists(os.path.join(root, f)):
        problem("missing figure:", f)
if "acmart" not in texts["paper.tex"].split("\n")[0]:
    problem("paper.tex does not use the acmart class")

# page limit of the built PDF: 12 pages excluding references
pdf = os.path.join(root, "paper.pdf")
if os.path.exists(pdf) and shutil.which("pdftotext"):
    pages = []
    i = 1
    while True:
        r = subprocess.run(["pdftotext", "-f", str(i), "-l", str(i), pdf, "-"], capture_output=True, text=True)
        if r.returncode != 0 or not r.stdout.strip():
            break
        pages.append(r.stdout)
        i += 1
    ref_page = next((k + 1 for k, pg in enumerate(pages) if re.search(r"^REFERENCES\s*$", pg, re.M | re.I)), None)
    print(f"paper.pdf: {len(pages)} pages, references start on page {ref_page}")
    body_ok = ref_page is not None and ref_page <= 12
    if ref_page == 13:   # allowed if page 13 holds nothing of the body before the reference list
        before = re.split(r"^REFERENCES\s*$", pages[12], maxsplit=1, flags=re.M | re.I)[0]
        lines = [ln.strip() for ln in before.splitlines()]
        lines = [ln for ln in lines if ln and not ln.isdigit() and not ln.startswith("Counters for Marginals")]
        body_ok = not lines
        print("page 13 before the references:", lines[:3] if lines else "nothing")
    if not body_ok:
        problem("the body of paper.pdf exceeds 12 pages (or the reference list was not found)")
print(f"{len(used)} number macros used, {sum(map(len, labels.values()))} labels, {len(cited)} citations; problems: {bad}")
sys.exit(1 if bad else 0)
