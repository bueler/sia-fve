
# for maint branch:
#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules

# for next branch:
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

mahaffy: mahaffy.o exactsia.o io.o q1op.o sia.o chkopts
	-${CLINKER} -o mahaffy mahaffy.o exactsia.o io.o q1op.o sia.o  ${PETSC_SNES_LIB}
	${RM} mahaffy.o exactsia.o io.o q1op.o sia.o

# for generating dome verification convergence figure in ../paper/
dome.verif: mahaffy domeconv.sh
	./domeconv.sh dome.verif 6 6
true.verif: mahaffy domeconv.sh
	./domeconv.sh true.verif 4 6 "-mah_true -snes_fd_color -da_grid_x 20 -da_grid_y 20"

# for generating bedstep verification table in ../paper/
bedstep.verif: mahaffy bedstepconv.sh
	./bedstepconv.sh bedstep.verif 4 6 "-snes_max_it 200 -mah_Neps 12"
bedstepnoup.verif: mahaffy bedstepconv.sh
	./bedstepconv.sh bedstepnoup.verif 4 6 "-snes_max_it 200 -mah_Neps 12 -mah_lambda 0.0"
bedstepfullup.verif: mahaffy bedstepconv.sh
	./bedstepconv.sh bedstepfullup.verif 4 6 "-snes_max_it 200 -mah_Neps 12 -mah_lambda 1.0"

# these targets are just for regression testing
regress/dome.short: mahaffy domeconv.sh
	./domeconv.sh regress/dome.short 1 6
regress/bedstep.short: mahaffy bedstepconv.sh
	./bedstepconv.sh regress/bedstep.short 0 6
regress/true.short: mahaffy domeconv.sh
	./domeconv.sh regress/true.short 1 6 "-mah_true -snes_fd_color -da_grid_x 15 -da_grid_y 15"
regress/noupwind.short: mahaffy domeconv.sh
	./domeconv.sh regress/noupwind.short 0 6 "-mah_lambda 0.0"
test: regress/dome.short regress/bedstep.short regress/true.short regress/noupwind.short
	git diff regress/

.PHONY: distclean texclean

dirclean:
	@rm -f *~ mahaffy *.pyc

