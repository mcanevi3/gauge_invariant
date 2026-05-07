@echo off
pdflatex root
bibtex root
pdflatex root
pdflatex root
del *.log
del *.aux
del *.out
