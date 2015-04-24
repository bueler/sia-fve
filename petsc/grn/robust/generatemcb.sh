#!/bin/bash

# this script generates a sequence of files with refining grids
#   mcb0.nc, mcb1.nc, ...
# based on the mass-conserving bed data

# do ./study.sh after this

if [ -z ${1+x} ]; then
  echo "ERROR: run as './generatemcb.sh N' with N=0,...,6"
  exit
fi

set -e # exit on error
set -x

if [ -z ${PISM+x} ]; then
  PISM=~/pism
fi

# MCB needs tools: NCO, ../inplace.py, and PISM util script nc2cdo.py
ncks -O -v x1,y1,thk,topg,climatic_mass_balance,mapping ../pism_Greenland_5km_v1.1.nc searise5km.nc
../inplace.py --fixcmbunits --ranges searise5km.nc
../inplace.py --oceanfix --ranges searise5km.nc
$PISM/util/nc2cdo.py searise5km.nc

# correspond to 9000m, 4500m, 3000m, 1500m, 1200m, 900m, 600m grids
BLOCKLEVELS="60 30 20 10 8 6 4"
COUNTER=0
for LEV in $BLOCKLEVELS; do
  ../mcb/average.py --block $LEV tmp${COUNTER}.nc
  ../mcb/remap2mcb.py searise5km.nc tmp${COUNTER}.nc mcb${COUNTER}.nc
  rm tmp${COUNTER}.nc
  if [ $COUNTER -eq $1 ]; then
    exit
  fi
  let COUNTER=COUNTER+1
done

set +x

echo ""
echo "********  DONE.  NOW RUN study.sh.  ********"

