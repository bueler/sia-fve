#!/bin/bash

set -e # exit on error
set -x

if [ -z ${1+x} ]; then
  NN=6
else
  NN=$1
fi

# move old result out of way if present
OUT=study.1500m
touch $OUT
mv -f $OUT SAFE_$OUT

# perform robustness study
STDOPTS="-cs_D0 30.0 -snes_monitor -snes_max_it 200 -pc_type asm -sub_pc_type lu"

NAMEROOT="mcb1500mS"

DAT=${NAMEROOT}4.dat
../../nc2petsc.py ${NAMEROOT}4.nc $DAT
DUMP=r1500m4/
mkdir -p $DUMP
mpiexec -n $NN ../../mahaffy -mah_read $DAT $STDOPTS -mah_notry -mah_dump $DUMP
head -n 1 ${DUMP}history.txt >> $OUT
cat ${DUMP}history.txt | grep "last successful value of eps" | sed 's/.* //g' >> $OUT
cat ${DUMP}history.txt | grep "maximum solution diffusivity" | sed 's/.* //g' >> $OUT
cat ${DUMP}history.txt | grep "total time" | sed 's/.* //g' >> $OUT

# runs below are time-stepping because
# we already know that continuation scheme fails on the steady state problem
# on the resolved-bed cases
LEVELS="3 2 1"

for LEV in $LEVELS; do
    DAT=${NAMEROOT}${LEV}.dat
    ../../nc2petsc.py ${NAMEROOT}${LEV}.nc $DAT
    PREV="-mah_readinitialsurface ${DUMP}unnamed.dat"
    DUMP=r1500m${LEV}/
    mkdir -p $DUMP
    mpiexec -n $NN ../../mahaffy -mah_T 1.0 -mah_dt 0.05 -cs_start 11 -mah_read $DAT $PREV $STDOPTS -mah_dump $DUMP -mah_showdata -draw_pause 5
    head -n 1 ${DUMP}history.txt >> $OUT
    cat ${DUMP}history.txt | grep "last successful value of eps" | sed 's/.* //g' >> $OUT
    cat ${DUMP}history.txt | grep "maximum solution diffusivity" | sed 's/.* //g' >> $OUT
    cat ${DUMP}history.txt | grep "total time" | sed 's/.* //g' >> $OUT
done

set +x

echo ""
echo "********  DONE.  SEE RESULT IN $OUT.  ********"

