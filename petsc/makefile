all: mahaffy

# for master branch:
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

base := appctx.o continuationscheme.o sia.o exactsia.o io.o q1op.o solver.o
base := $(addprefix base/, $(base))

mahaffy: mahaffy.o $(base) chkopts
	-${CLINKER} -o mahaffy mahaffy.o $(base) ${PETSC_LIB}
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
regress/unequal.short: mahaffy
	./mahaffy -da_grid_x 15 -da_grid_y 30 -mah_history tmp
	rm -f regress/unequal.short
	cat tmphistory.txt | grep "spacing in x"        | sed 's/.* //g' >> regress/unequal.short
	cat tmphistory.txt | grep "spacing in y"        | sed 's/.* //g' >> regress/unequal.short
	cat tmphistory.txt | grep "last successful value of eps" | sed 's/.* //g' >> regress/unequal.short
	cat tmphistory.txt | grep "max thickness error" | sed 's/.* //g' >> regress/unequal.short
	cat tmphistory.txt | grep "av thickness error"  | sed 's/.* //g' >> regress/unequal.short
	rm -f tmphistory.txt

test: regress/dome.short regress/bedstep.short regress/true.short regress/noupwind.short regress/unequal.short
	git diff regress/

.PHONY: dirclean

dirclean:
	@rm -f *~ base/*~ grn/*~ search/*~ regress/*~ verif/*~ grn/mcb/*~ mahaffy *.pyc grn/*.pyc search/*.pyc

