#!/bin/bash

set -e # exit on error
set -x

if [ -z ${1+x} ]; then
  NN=6
else
  NN=$1
fi

# move old result out of way if present
OUT=study.steadyfinemcb
touch $OUT
mv -f $OUT SAFE_$OUT

# perform robustness study
STDOPTS="-cs_D0 10.0 -snes_monitor -snes_max_it 400 -pc_type asm -sub_pc_type lu"

NAMEROOT="mcb1500mS"

DAT=${NAMEROOT}3.dat
../../nc2petsc.py ${NAMEROOT}3.nc $DAT
DUMP=result3/
mpiexec -n $NN ../../mahaffy -mah_read $DAT $STDOPTS -mah_notry -mah_dump $DUMP
cat $DAT >> $OUT
cat $DUMP/history.txt | grep "last successful value of eps" | sed 's/.* //g' >> $OUT
cat $DUMP/history.txt | grep "maximum solution diffusivity" | sed 's/.* //g' >> $OUT
cat $DUMP/history.txt | grep "total time" | sed 's/.* //g' >> $OUT

# runs below are time-stepping because
# we already know that continuation scheme fails on the steady state problem
# on the resolved-bed cases
LEVELS="2 1"
PREV=
for LEV in $LEVELS; do
    DAT=${NAMEROOT}${LEV}.dat
    ../../nc2petsc.py ${NAMEROOT}${LEV}.nc $DAT
    DUMP=result${LEV}/
    mpiexec -n $NN ../../mahaffy -cs_start 8 -mah_steps 10 -mah_dt 100.0 -mah_read $DAT $PREV $STDOPTS -mah_notry -mah_dump $DUMP
    cat $DAT >> $OUT
    cat $DUMP/history.txt | grep "last successful value of eps" | sed 's/.* //g' >> $OUT
    cat $DUMP/history.txt | grep "maximum solution diffusivity" | sed 's/.* //g' >> $OUT
    cat $DUMP/history.txt | grep "total time" | sed 's/.* //g' >> $OUT
    PREV="-mah_readinitial ${DUMP}/unnamed.dat"
done

set +x

echo ""
echo "********  DONE.  SEE RESULT IN $OUT.  ********"

