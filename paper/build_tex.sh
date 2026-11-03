#!/bin/bash
# Build the LaTeX source and PDF of the note with pandoc + pdflatex in the pandoc/latex Docker image (nothing installed
# on the host). Run from the repository root in Git Bash:  bash paper/build_tex.sh
# Outputs: paper/exact_graphstate_sampling_tex.md (converted Markdown), _tex.tex (the source arXiv compiles), _tex.pdf.
set -e
export MSYS_NO_PATHCONV=1
PYTHONIOENCODING=utf-8 python paper/md2tex.py paper/exact_graphstate_sampling.md
OPTS="-V geometry:margin=2.2cm -V fontsize=10pt -V colorlinks=true -V linkcolor=blue -V urlcolor=blue -V citecolor=blue"
docker run --rm -v "D:/GIT/Innovation/paper:/data" pandoc/latex:3.6 exact_graphstate_sampling_tex.md -s -o exact_graphstate_sampling_tex.tex $OPTS
docker run --rm -v "D:/GIT/Innovation/paper:/data" pandoc/latex:3.6 exact_graphstate_sampling_tex.md -o exact_graphstate_sampling_tex.pdf --pdf-engine=pdflatex $OPTS
ls -l paper/exact_graphstate_sampling_tex.tex paper/exact_graphstate_sampling_tex.pdf
