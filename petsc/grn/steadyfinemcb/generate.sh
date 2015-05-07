#!/bin/bash

# this script generates a sequence of files with refining grids
#   mcb0.nc,  mcb1.nc,  mcb2.nc,  mcb3.nc,  mcb4.nc
# at resolutions
#   9000m,    5400m,    3600m,    1800m,    900m
# based on the mass-conserving bed data, with blocking of
#   60,       36,       24,       12,       6
# from the 150m data

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

BLOCKLEVELS="60 36 24 12 6"
COUNTER=0
for LEV in $BLOCKLEVELS; do
  ../mcb/average.py --block $LEV tmp${COUNTER}.nc
  ../mcb/remap2mcb.py searise5km.nc tmp${COUNTER}.nc mcb${COUNTER}.nc
  rm tmp${COUNTER}.nc
  let COUNTER=COUNTER+1
done

set +x

echo ""
echo "********  DONE  ********"

