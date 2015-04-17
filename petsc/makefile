
# for maint branch:
#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules

# for master branch:
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

base := exactsia.o io.o q1op.o sia.o continuationscheme.o
base := $(addprefix base/, $(base))

mahaffy: mahaffy.o $(base) chkopts
	-${CLINKER} -o mahaffy mahaffy.o $(base) ${PETSC_SNES_LIB}
	${RM} mahaffy.o $(base)

# these targets are for regression testing
regress/dome.short: mahaffy verif/domeconv.sh
	verif/domeconv.sh regress/dome.short 1 6
regress/bedstep.short: mahaffy verif/bedstepconv.sh
	verif/bedstepconv.sh regress/bedstep.short 0 6
regress/true.short: mahaffy verif/domeconv.sh
	verif/domeconv.sh regress/true.short 1 6 "-mah_true -snes_fd_color -da_grid_x 15 -da_grid_y 15"
regress/noupwind.short: mahaffy verif/domeconv.sh
	verif/domeconv.sh regress/noupwind.short 0 6 "-mah_lambda 0.0"

test: regress/dome.short regress/bedstep.short regress/true.short regress/noupwind.short
	git diff regress/

.PHONY: distclean texclean

dirclean:
	@rm -f *~ mahaffy *.pyc

