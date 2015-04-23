#!/bin/bash

# this script generates two sequence of files
#   sea0.nc, sea1.nc, ...
#   fix0.nc, fix1.nc, ...
# for reading by ../../mahaffy

# do ../quickstart.sh first

# do ./study.sh after this

set -e # exit on error
set -x

if [ -z ${LONG+x} ]; then
  REFINELEVELS="2"
else
  # correspond to 2500m, 1667m, 1250m, 1000m, 625m grids, where refinement needed
  REFINELEVELS="2 3 4 5 8"
fi

# FIXME: check if ../grn exists

# SEARISE
# generate .nc
# 10 km:
ncks -v topg,cmb,thk -d x1,,,2 -d y1,,,2 ../grn.nc sea0.nc
# 5 km:
ncks -v topg,cmb,thk ../grn.nc sea1.nc
# 2.5 km and finer:
COUNTER=2
for LEV in $REFINELEVELS; do
  ../refine.py --factor $LEV ../grn.nc sea${COUNTER}.nc
  let COUNTER=COUNTER+1
done

if [ -z ${LONG+x} ]; then
  BLOCKLEVELS="60 30 20"
else
  # correspond to 9000m, 4500m, 3000m, 1500m, 1200m, 900m, 450m grids
  BLOCKLEVELS="60 30 20 10 8 6 3"
fi

# MCB
ncks -O -v x1,y1,thk,topg,climatic_mass_balance,mapping ../pism_Greenland_5km_v1.1.nc searise5km.nc
../inplace.py --fixcmbunits --ranges searise5km.nc
../inplace.py --oceanfix --ranges searise5km.nc
~/pism/util/nc2cdo.py searise5km.nc

COUNTER=0
for LEV in $BLOCKLEVELS; do
  ../mcb/average.py --block $LEV tmp${COUNTER}.nc
  ../mcb/remap2mcb.py searise5km.nc tmp${COUNTER}.nc mcb${COUNTER}.nc
  rm tmp${COUNTER}.nc
  let COUNTER=COUNTER+1
done

set +x

echo ""
echo "********  DONE.  NOW RUN study.sh.  ********"

