
# for maint branch:
#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules

# for next branch:
include ${PETSC_DIR}/lib/petsc-conf/variables
include ${PETSC_DIR}/lib/petsc-conf/rules

mahaffy: mahaffy.o exactsia.o binaryio.o chkopts
	-${CLINKER} -o mahaffy mahaffy.o exactsia.o binaryio.o  ${PETSC_SNES_LIB}
	${RM} mahaffy.o exactsia.o binaryio.o

.PHONY: distclean texclean

dirclean:
	@rm -f *~ mahaffy *.pyc

