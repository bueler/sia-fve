all: mahaffyfem.pdf

figures =
tikzfigs = fig1a.tikz fig1b.tikz fig2a.tikz fig2b.tikz fig3.tikz

mahaffyfem.pdf: mahaffyfem.aux mahaffyfem.tex $(figures) $(tikzfigs)
	pdflatex mahaffyfem
	pdflatex mahaffyfem

mahaffyfem.aux: mahaffyfem.tex $(figures) $(tikzfigs)
	pdflatex mahaffyfem
	bibtex mahaffyfem

mahaffyfem.bbl: mahaffyfem.aux ../paper/lc.bib
	bibtex mahaffyfem

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *~
