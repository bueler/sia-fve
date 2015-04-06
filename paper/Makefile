all: siafve.pdf

figures = domeverif.pdf bedstepverif.pdf domeprofile.pdf bedstepprofiles.pdf
tikzfigs = fig1a.tikz fig1b.tikz fig2a.tikz fig2b.tikz fig3.tikz

siafve.pdf: siafve.aux siafve.tex bedsteptable.tex $(figures) $(tikzfigs)
	pdflatex siafve
	pdflatex siafve

siafve.aux: siafve.tex bedsteptable.tex $(figures) $(tikzfigs)
	pdflatex siafve
	bibtex siafve

siafve.bbl: siafve.aux ../paper/lc.bib
	bibtex siafve

domeverif.pdf: convfigs.py ../petsc/dome.verif ../petsc/true.verif
	./convfigs.py --showdata -dome ../petsc/dome.verif -true ../petsc/true.verif -o domeverif.pdf

bedstepverif.pdf: convfigs.py ../petsc/bedstep.verif
	./convfigs.py --showdata -bedstep ../petsc/bedstep.verif -o bedstepverif.pdf

bedsteptable.tex: convfigs.py ../petsc/bedstep.verif
	./convfigs.py --bedsteptable -bedstep ../petsc/bedstep.verif > bedsteptable.tex

domeprofile.pdf:
	rm -rf tmpdome/
	mkdir tmpdome/
	../petsc/mahaffy -mah_dump tmpdome/ -da_refine 2 -snes_type vinewtonssls
	(cd tmpdome/ && ../../petsc/figsmahaffy.py --profile --half)
	cp tmpdome/HHexact1d.pdf domeprofile.pdf
	rm -rf tmpdome/

bedstepprofiles.pdf:
	rm -rf ../petsc/tmpbs/
	mkdir ../petsc/tmpbs/
	(cd ../petsc/ && ./bedstepprofiles.sh tmpbs/)
	cp ../petsc/tmpbs/HHexact1d.pdf bedstepprofiles.pdf
	rm -rf ../petsc/tmpbs/

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *~
