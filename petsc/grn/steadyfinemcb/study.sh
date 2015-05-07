#!/bin/bash

set -e # exit on error
set -x

if [ -z ${1+x} ]; then
  NN=6
else
  NN=$1
fi

if [ -z ${2+x} ]; then
  RES=1800
  START=4
  LEVELS="3 2 1"
else  # note argument 2 is *any* nonempty string
  RES=900
  START=5
  LEVELS="4 3 2 1"
fi

# move old result out of way if present
OUT=study.${RES}m
touch $OUT
mv -f $OUT SAFE_$OUT

cathistory ()
{
head -n 1 $1 >> $2
cat $1 | grep "last successful value of eps" | sed 's/.* //g' >> $2
cat $1 | grep "maximum solution diffusivity" | sed 's/.* //g' >> $2
cat $1 | grep "total time" | sed 's/.* //g' >> $2
}

# perform robustness study
STDOPTS="-cs_D0 30.0 -snes_monitor -snes_max_it 200 -pc_type asm -sub_pc_type lu"

NAMEROOT="mcb${RES}mS"

DAT=${NAMEROOT}${START}.dat
../../nc2petsc.py ${NAMEROOT}${START}.nc $DAT
DUMP=r${RES}m${START}/
mkdir -p $DUMP
mpiexec -n $NN ../../mahaffy -mah_read $DAT $STDOPTS -mah_notry -mah_dump $DUMP
cathistory ${DUMP}history.txt $OUT

# runs below are time-stepping because
# we already know that continuation scheme fails on the steady state problem
# on the resolved-bed cases
for LEV in $LEVELS; do
    DAT=${NAMEROOT}${LEV}.dat
    ../../nc2petsc.py ${NAMEROOT}${LEV}.nc $DAT
    PREV="-mah_readinitialsurface ${DUMP}unnamed.dat"
    DUMP=r${RES}m${LEV}/
    mkdir -p $DUMP
    mpiexec -n $NN ../../mahaffy -mah_T 1.0 -mah_dt 0.05 -cs_start 11 -mah_read $DAT $PREV $STDOPTS -mah_dump $DUMP
    # add if interactive: -mah_showdata -draw_pause 5
    cathistory ${DUMP}history.txt $OUT
done

set +x

echo ""
echo "********  DONE.  SEE RESULT IN $OUT.  ********"

