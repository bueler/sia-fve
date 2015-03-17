
# for maint branch:
#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules

# for next branch:
include ${PETSC_DIR}/lib/petsc-conf/variables
include ${PETSC_DIR}/lib/petsc-conf/rules

mahaffy: mahaffy.o exactsia.o io.o chkopts
	-${CLINKER} -o mahaffy mahaffy.o exactsia.o io.o  ${PETSC_SNES_LIB}
	${RM} mahaffy.o exactsia.o io.o

# these targets are used for generating verification figures in ../../mahaffypaper/
dome.verif: mahaffy domeconv.sh
	./domeconv.sh dome.verif 6 6
bedstep.verif: mahaffy bedstepconv.sh
	./bedstepconv.sh bedstep.verif 4 6

# these targets are just for regression testing
dome.short: mahaffy domeconv.sh
	./domeconv.sh dome.short 2 6
bedstep.short: mahaffy bedstepconv.sh
	./bedstepconv.sh bedstep.short 1 6

.PHONY: distclean texclean

dirclean:
	@rm -f *~ mahaffy *.pyc

