#!/usr/bin/env bash
# Build paper/paper.pdf and paper/appendix.pdf (acmart, sigconf). The two documents reference each
# other through the xr package, so they are compiled alternately until the references settle.
# `build_paper.sh camera` builds the camera-ready variant paper-camera.pdf instead (authors.tex, acks.tex;
# the submission files are not touched).
set -eu
cd "$(dirname "$0")/../paper"
if [ "${1:-}" = camera ]; then
  run() { pdflatex -interaction=nonstopmode -halt-on-error -jobname=paper-camera '\def\cameraready{}\input{paper}' > paper-camera.build.log 2>&1 || { tail -40 paper-camera.build.log; exit 1; }; }
  [ -s appendix.aux ] || { echo "build the submission first (appendix.aux is needed for the cross-references)"; exit 1; }
  run; bibtex paper-camera > /dev/null; run; run
  echo "paper-camera.pdf: $(pdfinfo paper-camera.pdf | awk '/^Pages/ {print $2}') pages"
  grep -E "Overfull \\\\hbox|undefined|Undefined" paper-camera.log | sort | uniq -c | head || true
  exit 0
fi
run() { pdflatex -interaction=nonstopmode -halt-on-error "$1" > "$1.build.log" 2>&1 || { tail -40 "$1.build.log"; exit 1; }; }
for doc in appendix paper; do run $doc; bibtex $doc > /dev/null; done
for i in 1 2; do run appendix; run paper; done
for doc in paper appendix; do
  pages=$(pdfinfo $doc.pdf 2>/dev/null | awk '/^Pages/ {print $2}')
  echo "$doc.pdf: ${pages:-?} pages"
  grep -E "Overfull \\\\hbox|undefined|Undefined|multiply defined" $doc.log | sort | uniq -c | head -20 || true
done
