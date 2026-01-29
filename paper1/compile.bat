@echo off
pdflatex autosam
pdflatex autosam
del *.log;
del *.aux;
del *.out;