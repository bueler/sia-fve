
# for maint branch:
#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules

# for next branch:
include ${PETSC_DIR}/lib/petsc-conf/variables
include ${PETSC_DIR}/lib/petsc-conf/rules

mahaffy: mahaffy.o exactsia.o io.o chkopts
	-${CLINKER} -o mahaffy mahaffy.o exactsia.o io.o  ${PETSC_SNES_LIB}
	${RM} mahaffy.o exactsia.o io.o

dome.verif: mahaffy domeconv.sh
	./domeconv.sh dome.verif 6 6

bedstep.verif: mahaffy bedstepconv.sh
	./bedstepconv.sh bedstep.verif 4 1

.PHONY: distclean texclean

dirclean:
	@rm -f *~ mahaffy *.pyc

