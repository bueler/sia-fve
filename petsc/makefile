
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

layer: layer.o  chkopts
	-${CLINKER} -o layer layer.o  ${PETSC_SNES_LIB}
	${RM} layer.o

mahaffy: mahaffy.o  chkopts
	-${CLINKER} -o mahaffy mahaffy.o  ${PETSC_SNES_LIB}
	${RM} mahaffy.o

doclayer.pdf: doclayer.tex
	pdflatex doclayer.tex
	pdflatex doclayer.tex

.PHONY: distclean texclean

dirclean:
	@rm -f *~ layer mahaffy *.pyc

texclean:
	@rm -f *~ *.aux *.log *.synctex.gz
