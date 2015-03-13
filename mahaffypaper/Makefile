all: mahaffyfem.pdf

figures = domeverif.pdf bedstepverif.pdf #bedstepprofile.pdf
tikzfigs = fig1a.tikz fig1b.tikz fig2a.tikz fig2b.tikz fig3.tikz

mahaffyfem.pdf: mahaffyfem.aux mahaffyfem.tex $(figures) $(tikzfigs)
	pdflatex mahaffyfem
	pdflatex mahaffyfem

mahaffyfem.aux: mahaffyfem.tex $(figures) $(tikzfigs)
	pdflatex mahaffyfem
	bibtex mahaffyfem

mahaffyfem.bbl: mahaffyfem.aux ../paper/lc.bib
	bibtex mahaffyfem

domeverif.pdf: convfigs.py ../petsc/mahaffy/dome.verif
	./convfigs.py -o domeverif.pdf ../petsc/mahaffy/dome.verif

bedstepverif.pdf: convfigs.py ../petsc/mahaffy/bedstep.verif
	./convfigs.py --bedstep -o bedstepverif.pdf ../petsc/mahaffy/bedstep.verif

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *~
