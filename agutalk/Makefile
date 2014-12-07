all: talk.pdf

figures = cartoon-layer.pdf cartoon-wclimate.pdf cartoon-sensitive-one.pdf cartoon-sensitive-two.pdf cartoon-sensitive-three.pdf cartoon-sets.pdf

cartoon-layer.pdf cartoon-wclimate.pdf cartoon-sensitive-one.pdf cartoon-sensitive-two.pdf cartoon-sensitive-three.pdf cartoon-sets.pdf:
	python figs1d.py
	#for X in layer wclimate sensitive-one; do pdfcrop cartoon-$$X.pdf cartoon-$$X.pdf; done
	pdfcrop cartoon-sets.pdf cartoon-sets.pdf

talk.aux: talk.tex $(figures)
	pdflatex talk

talk.pdf: talk.aux $(figures)
	pdflatex talk

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *.nav *.snm *.toc *~ cartoon-*.pdf
