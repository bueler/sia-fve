
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

advectlayer: advectlayer.o  chkopts
	-${CLINKER} -o advectlayer advectlayer.o  ${PETSC_SNES_LIB}
	${RM} advectlayer.o

.PHONY: distclean

distclean:
	@rm -f *~ advectlayer
