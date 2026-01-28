@echo off
pdflatex main 
pdflatex main 
del *.aux
del *.log
del *.out
