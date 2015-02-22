all: mahaffyfem.pdf

figures =
tikzfigs =

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
