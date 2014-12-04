all: talk.pdf

figures = layer-cartoon.pdf

layer-cartoon.pdf:
	python figs1d.py
	pdfcrop layer-cartoon.pdf layer-cartoon.pdf

talk.aux: talk.tex $(figures)
	pdflatex talk

talk.pdf: talk.aux $(figures)
	pdflatex talk

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *.nav *.snm *.toc *~
