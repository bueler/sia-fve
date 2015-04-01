all: siafve.pdf

figures = domeverif.pdf bedstepverif.pdf domeprofile.pdf bedstepprofile.pdf
tikzfigs = fig1a.tikz fig1b.tikz fig2a.tikz fig2b.tikz fig3.tikz

siafve.pdf: siafve.aux siafve.tex bedsteptable.tex $(figures) $(tikzfigs)
	pdflatex siafve
	pdflatex siafve

siafve.aux: siafve.tex bedsteptable.tex $(figures) $(tikzfigs)
	pdflatex siafve
	bibtex siafve

siafve.bbl: siafve.aux ../paper/lc.bib
	bibtex siafve

domeverif.pdf: convfigs.py ../petsc/dome.verif
	./convfigs.py -o domeverif.pdf ../petsc/dome.verif

bedstepverif.pdf: convfigs.py ../petsc/bedstep.verif
	./convfigs.py --bedstep -o bedstepverif.pdf ../petsc/bedstep.verif

bedsteptable.tex: convfigs.py ../petsc/bedstep.verif
	./convfigs.py --bedsteptable ../petsc/bedstep.verif > bedsteptable.tex

domeprofile.pdf:
	rm -rf tmpdome/
	mkdir tmpdome/
	../petsc/mahaffy -mah_dump tmpdome/ -da_refine 2 -snes_type vinewtonssls
	(cd tmpdome/ && ../../petsc/figsmahaffy.py --profile --half)
	cp tmpdome/HHexact1d.pdf domeprofile.pdf
	rm -rf tmpdome/

bedstepprofile.pdf:
	rm -rf tmpbs/
	mkdir tmpbs/
	mpiexec -n 4 ../petsc/mahaffy -mah_dump tmpbs/ -mah_bedstep -mah_D0 0.001 -da_refine 2 -snes_max_it 400 -snes_type vinewtonrsls -pc_type asm -sub_pc_type lu
	(cd tmpbs/ && ../../petsc/figsmahaffy.py --profile --half --sharpbed)
	cp tmpbs/HHexact1d.pdf bedstepprofile.pdf
	rm -rf tmpbs/

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *~
