all: talk.pdf

figures = cartoon-layer.pdf cartoon-wclimate.pdf cartoon-sensitive-one.pdf cartoon-sensitive-two.pdf cartoon-sensitive-three.pdf cartoon-sets.pdf

cartoon-layer.pdf cartoon-wclimate.pdf cartoon-sensitive-one.pdf cartoon-sensitive-two.pdf cartoon-sensitive-three.pdf cartoon-sets.pdf: figs1d.py
	python figs1d.py
	pdfcrop cartoon-sets.pdf cartoon-sets.pdf

anim/u-1.png: genfigs.py ../petsc/layer.c
	(cd ../petsc/ && make layer)
	../petsc/layer -lay_genfigs anim/ -da_refine 3 -lay_dt 0.01 -lay_steps 100 -lay_adscheme 1
	python genfigs.py anim/ --stride 4

talk.aux: talk.tex $(figures)
	pdflatex talk

talk.pdf: talk.aux $(figures) anim/u-1.png
	pdflatex talk

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *.nav *.snm *.toc *~ cartoon-*.pdf anim/f.txt anim/u-*txt anim/u-*png
