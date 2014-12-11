
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

layer: layer.o  chkopts
	-${CLINKER} -o layer layer.o  ${PETSC_SNES_LIB}
	${RM} layer.o

.PHONY: distclean

distclean:
	@rm -f *~ layer
