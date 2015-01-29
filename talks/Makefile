all: agutalk.pdf liwgtalk.pdf

figures = cartoon-wclimate.pdf cartoon-sensitive-one.pdf cartoon-sensitive-two.pdf cartoon-sensitive-three.pdf cartoon-sets.pdf masstimeseries.pdf

cartoon-wclimate.pdf cartoon-sensitive-one.pdf cartoon-sensitive-two.pdf cartoon-sensitive-three.pdf cartoon-sets.pdf: figs1d.py
	python figs1d.py
	pdfcrop cartoon-sets.pdf cartoon-sets.pdf

../petsc/layer: ../petsc/layer.c
	(cd ../petsc/ && make layer)

anim/adv/u-1.png: genfigs.py ../petsc/layer
	(cd anim/ && rm -rf adv/ && mkdir adv/)
	../petsc/layer -lay_genfigs anim/adv/ -da_refine 3 -lay_dt 0.04 -lay_steps 25 -lay_adscheme 1 -lay_jac
	python genfigs.py anim/adv/ --stride 1 --uscale 12.0

anim/sia/u-1.png: genfigs.py ../petsc/layer
	(cd anim/ && rm -rf sia/ && mkdir sia/)
	../petsc/layer -lay_genfigs anim/sia/ -da_refine 3 -lay_dt 0.02 -lay_steps 500 -lay_adscheme 1 -lay_lambda 1.0 -lay_gamma 0.1 -lay_glenn 3.0 -lay_jac
	python genfigs.py anim/sia/ --stride 20 --uscale 0.9

mass0.txt:
	./refinemass.sh

masstimeseries.pdf: mass0.txt
	python ./masstimefig.py -o masstimeseries.pdf mass?.txt

agutalk.aux: agutalk.tex $(figures) anim/adv/u-1.png anim/sia/u-1.png
	pdflatex agutalk

agutalk.pdf: agutalk.aux $(figures) anim/adv/u-1.png anim/sia/u-1.png
	pdflatex agutalk


liwgtalk.aux: liwgtalk.tex $(figures) anim/adv/u-1.png anim/sia/u-1.png
	pdflatex liwgtalk

liwgtalk.pdf: liwgtalk.aux $(figures) anim/adv/u-1.png anim/sia/u-1.png
	pdflatex liwgtalk


.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *.nav *.snm *.toc *~ cartoon-*.pdf anim/f.txt anim/u-*txt anim/u-*png
