#!/bin/bash

# this script generates a sequence of files with refining grids
#   mcb1.nc, mcb2.nc, mcb3.nc, mcb4.nc
# at resolutions
#   4500m,   3000m,   1500m,   750m
# based on the mass-conserving bed data

# do ./setup1500m.sh or ./setup750m.sh after this

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

# correspond to 4500m, 3000m, 1500m, 750m grids
BLOCKLEVELS="30 20 10 5"
COUNTER=1
for LEV in $BLOCKLEVELS; do
  ../mcb/average.py --block $LEV tmp${COUNTER}.nc
  ../mcb/remap2mcb.py searise5km.nc tmp${COUNTER}.nc mcb${COUNTER}.nc
  rm tmp${COUNTER}.nc
  let COUNTER=COUNTER+1
done

set +x

echo ""
echo "********  DONE  ********"

