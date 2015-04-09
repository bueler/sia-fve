all: siafve.pdf

figures = domeverif.pdf bedstepverif.pdf domeprofile.pdf bedstepprofiles.pdf newtonconv.pdf
tikzfigs = fig1a.tikz fig1b.tikz fig2a.tikz fig2b.tikz fig3.tikz

siafve.pdf: siafve.bbl siafve.aux siafve.tex bedsteptable.tex $(figures) $(tikzfigs)
	pdflatex siafve
	pdflatex siafve

siafve.aux: siafve.tex bedsteptable.tex $(figures) $(tikzfigs)
	pdflatex siafve
	bibtex siafve

siafve.bbl: siafve.bib
	bibtex siafve

domeverif.pdf: convfigs.py ../petsc/dome.verif ../petsc/true.verif
	./convfigs.py --showdata -dome ../petsc/dome.verif -true ../petsc/true.verif -o domeverif.pdf

bedstepverif.pdf: convfigs.py ../petsc/bedstep.verif
	./convfigs.py --showdata -bedstep ../petsc/bedstep.verif -o bedstepverif.pdf

bedsteptable.tex: convfigs.py ../petsc/bedstep.verif ../petsc/bedstepnoup.verif ../petsc/bedstepfullup.verif
	./convfigs.py --bedsteptable -tablefiles ../petsc/bedstep.verif,../petsc/bedstepnoup.verif,../petsc/bedstepfullup.verif > bedsteptable.tex

newtonconv.pdf: newtonconv.txt
	./newtonconvshow.py -o newtonconv.pdf newtonconv.txt

newtonconv.txt:
	../petsc/mahaffy -da_refine 2 -mah_D0 1.0 -snes_type vinewtonssls -snes_monitor -snes_rtol 1.0e-10 > newtonconv.txt

domeprofile.pdf:
	rm -rf ../petsc/tmpdome/
	mkdir ../petsc/tmpdome/
	(cd ../petsc/ && ./domeprofileinset.sh tmpdome/)
	cp ../petsc/tmpdome/HHexact1d.pdf domeprofile.pdf
	rm -rf ../petsc/tmpdome/

bedstepprofiles.pdf:
	rm -rf ../petsc/tmpbs/
	mkdir ../petsc/tmpbs/
	(cd ../petsc/ && ./bedstepprofiles.sh tmpbs/)
	cp ../petsc/tmpbs/HHexact1d.pdf bedstepprofiles.pdf
	rm -rf ../petsc/tmpbs/

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *~
