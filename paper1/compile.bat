@echo off
pdflatex autosam
pdflatex autosam
bibtex autosam
pdflatex autosam
del *.log;
del *.aux;
del *.out;