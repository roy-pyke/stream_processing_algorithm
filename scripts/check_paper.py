#!/usr/bin/env python3
"""Static checks of the LaTeX sources (no TeX installation needed).

  - every \\NumXxx macro used in paper/*.tex is defined in paper/numbers.tex and is not "??"
  - every \\ref / \\cite target exists
  - \\begin / \\end environments and braces are balanced
  - every \\input and \\includegraphics file exists
Exit status 1 if anything is wrong.
"""
import os
import re
import sys

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper")
bad = 0
texts = {}
for name in ("paper.tex", "results.tex"):
    texts[name] = open(os.path.join(root, name)).read()
defs = dict(re.findall(r"\\newcommand\{\\(Num[A-Za-z]+)\}\{(.*)\}", open(os.path.join(root, "numbers.tex")).read()))
allt = "\n".join(texts.values())
used = sorted(set(re.findall(r"\\(Num[A-Za-z]+)", allt)))
for u in used:
    if u not in defs:
        print("UNDEFINED macro:", u)
        bad += 1
    elif "??" in defs[u] or "nan" in defs[u].lower():
        print("macro without value:", u, "=", defs[u])
        bad += 1
labels = set(re.findall(r"\\label\{([^}]+)\}", allt))
for r in sorted(set(re.findall(r"\\(?:eq)?ref\{([^}]+)\}", allt))):
    if r not in labels:
        print("undefined reference:", r)
        bad += 1
bib = set(re.findall(r"\\bibitem\{([^}]+)\}", allt))
cited = set()
for c in re.findall(r"\\cite\{([^}]+)\}", allt):
    cited.update(x.strip() for x in c.split(","))
for c in sorted(cited - bib):
    print("undefined citation:", c)
    bad += 1
for c in sorted(bib - cited):
    print("note: bibliography entry never cited:", c)
for name, t in texts.items():
    stack = []
    for m in re.finditer(r"\\(begin|end)\{([^}]+)\}", t):
        if m.group(1) == "begin":
            stack.append(m.group(2))
        elif not stack or stack.pop() != m.group(2):
            print(f"{name}: unbalanced environment near offset {m.start()}: {m.group(0)}")
            bad += 1
    if stack:
        print(f"{name}: unclosed environments: {stack}")
        bad += 1
    body = re.sub(r"(?<!\\)%.*", "", t)
    body = body.replace("\\{", "").replace("\\}", "")
    if body.count("{") != body.count("}"):
        print(f"{name}: unbalanced braces: {body.count('{')} open, {body.count('}')} close")
        bad += 1
    if len(re.findall(r"(?<!\\)\$", body)) % 2:
        print(f"{name}: odd number of $ signs")
        bad += 1
for f in re.findall(r"\\input\{([^}]+)\}", allt):
    path = os.path.join(root, f if f.endswith(".tex") else f + ".tex")
    if not os.path.exists(path):
        print("missing \\input file:", f)
        bad += 1
for f in re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", allt):
    if not os.path.exists(os.path.join(root, f)):
        print("missing figure:", f)
        bad += 1
print(f"{len(used)} number macros used, {len(labels)} labels, {len(cited)} citations; problems: {bad}")
sys.exit(1 if bad else 0)
