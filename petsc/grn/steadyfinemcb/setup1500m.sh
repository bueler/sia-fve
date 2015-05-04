#!/bin/bash

# generate 4 different versions of the 1500m MCB bed, with
#   S4 interpolated from the 9000m grid
#   S3 interpolated from the 4500m grid
#   S2 interpolated from the 3000m grid
#   S1 a copy of the 1500m grid

set -e # exit on error
set -x

../refine.py --factor 6 mcb0.nc mcb1500mS4.nc
../refine.py --factor 3 mcb1.nc mcb1500mS3.nc
../refine.py --factor 2 mcb2.nc mcb1500mS2.nc
../refine.py --factor 1 mcb3.nc mcb1500mS1.nc
for LEV in 1 2 3 4; do ncdump -h mcb1500mS${LEV}.nc |head -n 5; done

# next: chop so have same dims
for LEV in 1 2 3 4; do ncks -O -d x1,,990 -d y1,,1788 mcb1500mS${LEV}.nc mcb1500mS${LEV}.nc; done
for LEV in 1 2 3 4; do ncdump -h mcb1500mS${LEV}.nc |head -n 5; done

set +x

