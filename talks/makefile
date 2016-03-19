all: coppertalk.pdf

figures = cartoon-wclimate.pdf cartoon-sensitive-one.pdf cartoon-sensitive-two.pdf cartoon-sensitive-three.pdf cartoon-sets.pdf obstacle1d.pdf ncp1d.pdf vi1d.pdf rseps.pdf

cartoon-wclimate.pdf cartoon-sensitive-one.pdf cartoon-sensitive-two.pdf cartoon-sensitive-three.pdf cartoon-sets.pdf: figs1d.py
	python figs1d.py
	pdfcrop cartoon-sets.pdf cartoon-sets.pdf

obstacle1d.pdf ncp1d.pdf vi1d.pdf: obstacle1d.py
	python obstacle1d.py
	pdfcrop obstacle1d.pdf obstacle1d.pdf
	pdfcrop ncp1d.pdf ncp1d.pdf
	pdfcrop vi1d.pdf vi1d.pdf

rseps.pdf: rseps.py ../paper/study.robust
	./rseps.py -o rseps.pdf ../paper/study.robust

coppertalk.aux: coppertalk.tex $(figures)
	pdflatex coppertalk

coppertalk.pdf: coppertalk.aux $(figures)
	pdflatex coppertalk


.PHONY: clean

clean:
	@rm -f *~ *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz* *.nav *.snm *.toc *cartoon*.pdf *1d.pdf
