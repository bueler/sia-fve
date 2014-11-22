all: lc.pdf

figures = domains-fig.pdf

domains-fig.pdf:
	python domainsfigure.py
	pdfcrop domains-fig.pdf domains-fig.pdf

lc.pdf: lc.aux lc.bbl lc.tex $(figures)
	pdflatex lc

lc.aux: lc.tex $(figures)
	pdflatex lc
	bibtex lc

lc.bbl: lc.aux lc.bib
	bibtex lc

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *~
