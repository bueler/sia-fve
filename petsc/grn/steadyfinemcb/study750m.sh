#!/bin/bash

set -e # exit on error
set -x

if [ -z ${1+x} ]; then
  NN=6
else
  NN=$1
fi

# move old result out of way if present
OUT=study.750m
touch $OUT
mv -f $OUT SAFE_$OUT

# perform robustness study
STDOPTS="-cs_D0 30.0 -snes_monitor -snes_max_it 400 -pc_type asm -sub_pc_type lu"

NAMEROOT="mcb750mS"

DAT=${NAMEROOT}5.dat
../../nc2petsc.py ${NAMEROOT}5.nc $DAT
DUMP=r750m5/
mkdir -p $DUMP
mpiexec -n $NN ../../mahaffy -mah_read $DAT $STDOPTS -mah_notry -mah_dump $DUMP
head -n 1 $DUMP/history.txt >> $OUT
cat $DUMP/history.txt | grep "last successful value of eps" | sed 's/.* //g' >> $OUT
cat $DUMP/history.txt | grep "maximum solution diffusivity" | sed 's/.* //g' >> $OUT
cat $DUMP/history.txt | grep "total time" | sed 's/.* //g' >> $OUT

# runs below are time-stepping because we already know that continuation scheme
# fails on the steady state problem on the highly-resolved-bed
LEVELS="4 3 2 1"

for LEV in $LEVELS; do
    DAT=${NAMEROOT}${LEV}.dat
    ../../nc2petsc.py ${NAMEROOT}${LEV}.nc $DAT
    PREV="-mah_readinitialsurface ${DUMP}unnamed.dat"
    DUMP=r750m${LEV}/
    mkdir -p $DUMP
    # two-week time steps for 1.0 year
    mpiexec -n $NN ../../mahaffy -mah_T 1.0 -mah_dt 0.04166666 -cs_start 8 -mah_read $DAT $PREV $STDOPTS -mah_dump $DUMP
    head -n 1 $DUMP/history.txt >> $OUT
    cat $DUMP/history.txt | grep "last successful value of eps" | sed 's/.* //g' >> $OUT
    cat $DUMP/history.txt | grep "maximum solution diffusivity" | sed 's/.* //g' >> $OUT
    cat $DUMP/history.txt | grep "total time" | sed 's/.* //g' >> $OUT
done

set +x

echo ""
echo "********  DONE.  SEE RESULT IN $OUT.  ********"

