@echo off

pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex

del *.aux 
del *.log
del *.out
del *.toc