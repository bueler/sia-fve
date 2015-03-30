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

../petsc/mahaffy: ../petsc/mahaffy.c
	(cd ../petsc/ && make mahaffy)

domeprofile.pdf: ../petsc/mahaffy ../petsc/figsmahaffy.py
	rm -rf tmpdome/
	mkdir tmpdome/
	../petsc/mahaffy -mah_dump tmpdome/ -da_refine 2 -snes_type vinewtonssls
	(cd tmpdome/ && ../../petsc/figsmahaffy.py --profile --half)
	cp tmpdome/HHexact1d.pdf domeprofile.pdf
	rm -rf tmpdome/

bedstepprofile.pdf: ../petsc/mahaffy ../petsc/figsmahaffy.py
	rm -rf tmpbs/
	mkdir tmpbs/
	../petsc/mahaffy -mah_dump tmpbs/ -da_refine 2 -mah_bedstep -mah_D0 0.01
	(cd tmpbs/ && ../../petsc/figsmahaffy.py --profile --half --sharpbed)
	cp tmpbs/HHexact1d.pdf bedstepprofile.pdf
	rm -rf tmpbs/

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *~
