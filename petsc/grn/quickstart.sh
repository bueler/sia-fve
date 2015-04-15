#!/bin/bash


### RUN EXAMPLES ###
#  $ ./quickstart.sh
#  $ QSSHOW= QSREFINE=2 ./quickstart.sh "-mah_dtBE 10.0 -snes_max_it 200"

set -e # exit on error
set -x # echo commands

ncks -O -v x1,y1,thk,topg,climatic_mass_balance pism_Greenland_5km_v1.1.nc grn.nc

./inplace.py --fixcmbunits --ranges grn.nc

./inplace.py --oceanfix --ranges grn.nc

if [ -z ${QSMPI+x} ]; then
  QSMPI=6
fi

if [ -z ${QSREFINE+x} ]; then
  QSREFINE=1
fi

if [ -z ${QSSHOW+x} ]; then
  QSSHOW="-mah_showdata -draw_pause 2"
fi

if [ -z ${QSSOLVE+x} ]; then
  QSSOLVE="-pc_type asm -sub_pc_type lu"
fi

if [ -z ${QSD0+x} ]; then
  QSD0="10.0"
fi

if [ "$QSREFINE" -eq "1" ]; then
  NCDATA=grn.nc
  DATA=grn.dat
  DUMP=test/
else
  NCDATA=grnrefine$QSREFINE.nc
  DATA=grnrefine$QSREFINE.dat
  DUMP=testrefine$QSREFINE/
  ./refine.py --factor $QSREFINE grn.nc $NCDATA
fi

### OPTIONAL PROCESSING STEPS HERE ###
# coarsen:
#ncks -d x1,,,4 -d y1,,,4 grn.nc grn20km.nc
# apply bed sweeps:
#cp grn.nc grn_origtopg.nc
#./inplace.py --bedundip --sweeps 2 --ranges grn.nc
# trim so that -snes_fd_color works:
#ncks -O -d x1,0,299 grn.nc grn.nc

./nc2petsc.py $NCDATA $DATA

mkdir -p $DUMP

mpiexec -n $QSMPI ../mahaffy -mah_read $DATA -mah_D0 $QSD0 -snes_monitor $QSSHOW $QSSOLVE -mah_dump $DUMP $1

