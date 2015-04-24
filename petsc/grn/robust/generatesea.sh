#!/bin/bash

# this script generates a sequence of files with refining grids
#   sea0.nc, sea1.nc, ...
# based on the SeaRISE data

# do ../quickstart.sh first

# do ./study.sh after this

if [ -z ${1+x} ]; then
  echo "ERROR: run as './generatesea.sh N' with N=2,...,6"
  exit
fi

set -e # exit on error
set -x

# FIXME: check if ../grn exists

# 10 km:
ncks -v topg,cmb,thk -d x1,,,2 -d y1,,,2 ../grn.nc sea0.nc
# 5 km:
ncks -v topg,cmb,thk ../grn.nc sea1.nc
# correspond to 2500m, 1667m, 1250m, 1000m, 625m grids where refinement needed:
REFINELEVELS="2 3 4 5 8"
COUNTER=2
for LEV in $REFINELEVELS; do
  ../refine.py --factor $LEV ../grn.nc sea${COUNTER}.nc
  if [ $COUNTER -eq $1 ]; then
    exit
  fi
  let COUNTER=COUNTER+1
done

set +x

echo ""
echo "********  DONE.  ********"

