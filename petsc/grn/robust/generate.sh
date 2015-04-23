#!/bin/bash

# do ../quickstart.sh first

# FIXME: make optional depending on whether LONG is set
#REFINELEVELS="2 3 5 8"
REFINELEVELS="2"

set -e # exit on error
set -x

# FIXME: check if ../grn exists

# generate .nc
# 10 km:
ncks -d x1,,,2 -d y1,,,2 ../grn.nc sea0.nc
# 5 km:
cp ../grn.nc sea1.nc
# 2.5 km and finer:
for LEV in $REFINELEVELS; do
  ../refine.py --factor $LEV ../grn.nc sea${LEV}.nc
done

# FIXME: also generate mcb${LEV}.nc

set +x

echo ""
echo "********  DONE.  NOW RUN study.sh.  ********"

