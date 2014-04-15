all: lc.pdf

figures = 

# TO BUILD lc.pdf:
#   * ice-bib.bib needs to be a link to the same object in pism-dev/doc/
lc.pdf: lc.aux lc.bbl lc.tex $(figures)
	pdflatex lc

lc.aux: lc.tex $(figures)
	pdflatex lc
	bibtex lc

lc.bbl: lc.aux ice-bib.bib
	bibtex lc

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *~
