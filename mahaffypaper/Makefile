all: mahaffyfem.pdf

figures = domeverif.pdf bedstepverif.pdf domeprofile.pdf bedstepprofile.pdf
tikzfigs = fig1a.tikz fig1b.tikz fig2a.tikz fig2b.tikz fig3.tikz

mahaffyfem.pdf: mahaffyfem.aux mahaffyfem.tex bedsteptable.tex $(figures) $(tikzfigs)
	pdflatex mahaffyfem
	pdflatex mahaffyfem

mahaffyfem.aux: mahaffyfem.tex bedsteptable.tex $(figures) $(tikzfigs)
	pdflatex mahaffyfem
	bibtex mahaffyfem

mahaffyfem.bbl: mahaffyfem.aux ../paper/lc.bib
	bibtex mahaffyfem

domeverif.pdf: convfigs.py ../petsc/mahaffy/dome.verif
	./convfigs.py -o domeverif.pdf ../petsc/mahaffy/dome.verif

bedstepverif.pdf: convfigs.py ../petsc/mahaffy/bedstep.verif
	./convfigs.py --bedstep -o bedstepverif.pdf ../petsc/mahaffy/bedstep.verif

bedsteptable.tex: convfigs.py ../petsc/mahaffy/bedstep.verif
	./convfigs.py --bedsteptable ../petsc/mahaffy/bedstep.verif > bedsteptable.tex

../petsc/mahaffy/mahaffy: ../petsc/mahaffy/mahaffy.c
	(cd ../petsc/mahaffy/ && make mahaffy)

domeprofile.pdf: ../petsc/mahaffy/mahaffy ../petsc/mahaffy/figsmahaffy.py
	rm -rf tmpdome/
	mkdir tmpdome/
	../petsc/mahaffy/mahaffy -mah_dump tmpdome/ -da_refine 2
	(cd tmpdome/ && ../../petsc/mahaffy/figsmahaffy.py --profile --half)
	cp tmpdome/HHexact1d.pdf domeprofile.pdf
	rm -rf tmpdome/

bedstepprofile.pdf: ../petsc/mahaffy/mahaffy ../petsc/mahaffy/figsmahaffy.py
	rm -rf tmpbs/
	mkdir tmpbs/
	../petsc/mahaffy/mahaffy -mah_dump tmpbs/ -da_refine 2 -mah_bedstep -mah_D0 0.01
	(cd tmpbs/ && ../../petsc/mahaffy/figsmahaffy.py --profile --half --sharpbed)
	cp tmpbs/HHexact1d.pdf bedstepprofile.pdf
	rm -rf tmpbs/

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *~
