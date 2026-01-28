#!/bin/zsh
set -e

pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex

rm -f *.aux 
rm -f *.log
rm -f *.out
rm -f *.toc