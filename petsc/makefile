
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

layer: layer.o  chkopts
	-${CLINKER} -o layer layer.o  ${PETSC_SNES_LIB}
	${RM} layer.o

doclayer.pdf: doclayer.tex
	pdflatex doclayer.tex
	pdflatex doclayer.tex

.PHONY: distclean texclean

distclean:
	@rm -f *~ layer

texclean:
	@rm -f *~ *.aux *.log *.synctex.gz
