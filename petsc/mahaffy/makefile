
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
regress/dome.short: mahaffy domeconv.sh
	./domeconv.sh regress/dome.short 1 6
regress/bedstep.short: mahaffy bedstepconv.sh
	./bedstepconv.sh regress/bedstep.short 0 6
regress/true.short: mahaffy domeconv.sh
	./domeconv.sh regress/true.short 0 6 -mah_true
regress/noupwind.short: mahaffy domeconv.sh
	./domeconv.sh regress/noupwind.short 0 6 -mah_noupwind
test: regress/dome.short regress/bedstep.short regress/true.short regress/noupwind.short
	git diff regress/

.PHONY: distclean texclean

dirclean:
	@rm -f *~ mahaffy *.pyc

