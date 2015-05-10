#!/bin/bash

set -e # exit on error
set -x

if [ -z ${1+x} ]; then
  NN=6
else  # this case if a first argument is given
  NN=$1
fi

if [ -z ${FINE+x} ]; then
  RES=1800
  START=4
  LEVELS="3 2 1"
else  # this case if FINE is *any* nonempty string
  RES=900
  START=5
  LEVELS="4 3 2 1"
fi

if [ -z ${RUNT+x} ]; then
  RUNT=1.0
fi

if [ -z ${RUNDT+x} ]; then
  RUNDT=0.05
fi

EXEC=../../mahaffy
STDOPTS="-cs_D0 30.0 -snes_monitor -pc_type asm -sub_pc_type lu"
NAMEROOT="mcb${RES}mS"

# compute steady-state solution on highly-averaged bed
DAT=${NAMEROOT}${START}.dat
../../nc2petsc.py ${NAMEROOT}${START}.nc $DAT
DUMP=r${RES}m${START}/
mkdir -p $DUMP
mpiexec -n $NN $EXEC -mah_read $DAT $STDOPTS -snes_max_it 200 -mah_notry -mah_dump $DUMP

# runs below are time-stepping because we already know that continuation scheme
# fails on the steady state problem on the resolved-bed cases
for LEV in $LEVELS; do
    DAT=${NAMEROOT}${LEV}.dat
    ../../nc2petsc.py ${NAMEROOT}${LEV}.nc $DAT
    PREV=${DUMP}unnamed.dat
    DUMP=r${RES}m${LEV}/
    mkdir -p $DUMP
    mpiexec -n $NN $EXEC -mah_read $DAT -mah_readinitialsurface $PREV \
        $STDOPTS -mah_T $RUNT -mah_dt $RUNDT -cs_start 11 -mah_dump $DUMP
    # add if interactive: -mah_showdata -draw_pause 5
done

set +x

echo ""
echo "********  DONE.  ********"

