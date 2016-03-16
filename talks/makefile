all:

figures = cartoon-wclimate.pdf cartoon-sensitive-one.pdf cartoon-sensitive-two.pdf cartoon-sensitive-three.pdf cartoon-sets.pdf

# make link to layer-conserve/talks/figs1d.py
cartoon-wclimate.pdf cartoon-sensitive-one.pdf cartoon-sensitive-two.pdf cartoon-sensitive-three.pdf cartoon-sets.pdf: figs1d.py
	python figs1d.py
	pdfcrop cartoon-sets.pdf cartoon-sets.pdf

coppertalk.aux: coppertalk.tex $(figures)
	pdflatex coppertalk

coppertalk.pdf: coppertalk.aux $(figures)
	pdflatex coppertalk


.PHONY: clean

clean:
	@rm -f *~ *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz* *.nav *.snm *.toc *cartoon*.pdf
