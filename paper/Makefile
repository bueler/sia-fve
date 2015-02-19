all: lc.pdf

figures = domains-fig.pdf
tikzfigs = fvmesh.tikz

domains-fig.pdf: domainsfigure.py
	python domainsfigure.py
	pdfcrop domains-fig.pdf domains-fig.pdf

fvmesh.tikz: octagon.poly tri2tikz.py
	triangle -pqa0.085 octagon
	# FIXME this version is temporary: labels on nodes help development of removal alg
	python tri2tikz.py --nodesize 0.0 --scale 8.0 --noboundary --labelnodes octagon.1 fvmesh.tikz

lc.pdf: lc.aux lc.bbl lc.tex $(figures) $(tikzfigs)
	pdflatex lc
	pdflatex lc

lc.aux: lc.tex $(figures) $(tikzfigs)
	pdflatex lc
	bibtex lc

lc.bbl: lc.aux lc.bib
	bibtex lc

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *~
