all: siafve.pdf

figures = domeverif.pdf bedstepverif.pdf domeprofile.pdf bedstepprofiles.pdf newtonconv.pdf grnrobusteps.pdf
tikzfigs = fig1a.tikz fig1b.tikz fig2a.tikz fig2b.tikz fig3.tikz

siafve.pdf: siafve.aux
	bibtex siafve
	pdflatex siafve
	pdflatex siafve
	pdflatex siafve

siafve.aux: siafve.tex bedsteptable.tex $(figures) $(tikzfigs)
	pdflatex siafve

VERIF = ../petsc/verif/

domeverif.pdf: convfigs.py $(VERIF)dome.verif $(VERIF)true.verif
	./convfigs.py --showdata -dome $(VERIF)dome.verif -true $(VERIF)true.verif -o domeverif.pdf

bedstepverif.pdf: convfigs.py $(VERIF)bedstep.verif
	./convfigs.py --showdata -bedstep $(VERIF)bedstep.verif -o bedstepverif.pdf

bedsteptable.tex: convfigs.py $(VERIF)bedstep.verif $(VERIF)bedstepnoup.verif $(VERIF)bedstepfullup.verif
	./convfigs.py --bedsteptable -tablefiles $(VERIF)bedstep.verif,$(VERIF)bedstepnoup.verif,$(VERIF)bedstepfullup.verif > bedsteptable.tex

newtonconv.pdf: newtonconv.txt
	./newtonconvshow.py -o newtonconv.pdf newtonconv.txt

newtonconv.txt:
	../petsc/mahaffy -da_refine 2 -mah_D0 1.0 -snes_type vinewtonssls -snes_monitor -snes_rtol 1.0e-10 > newtonconv.txt

grnrobusteps.pdf: study.robust
	./robustfigs.py -o grnrobusteps.pdf study.robust

TMP = ../petsc/tmp/

domeprofile.pdf:
	rm -rf $(TMP)
	mkdir $(TMP)
	(cd ../petsc/ && ./domeprofileinset.sh tmp/)
	cp $(TMP)HHexact1d.pdf domeprofile.pdf
	rm -rf $(TMP)

bedstepprofiles.pdf:
	rm -rf $(TMP)
	mkdir $(TMP)
	(cd ../petsc/ && ./bedstepprofiles.sh tmp/)
	cp $(TMP)HHexact1d.pdf bedstepprofiles.pdf
	rm -rf $(TMP)

.PHONY: clean

clean:
	@rm -f *.pyc *.out *.aux *.log *.bbl *.blg *.synctex.gz *~

