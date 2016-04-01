include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

fractal: fractal.o  chkopts
	-${CLINKER} -o fractal fractal.o  ${PETSC_LIB}
	${RM} fractal.o

# etc

.PHONY: distclean

distclean:
	@rm -f *~ fractal

