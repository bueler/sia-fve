#!/bin/bash

# generate 4 different versions of the 750m MCB bed, with
#   S4 smoothed from the 4500m grid
#   S3 smoothed from the 3000m grid
#   S2 smoothed from the 1500m grid
#   S1 a copy of the 750m grid

set -e # exit on error
set -x

../refine.py --factor 6 mcb1.nc mcb750mS4.nc
../refine.py --factor 4 mcb2.nc mcb750mS3.nc
../refine.py --factor 2 mcb3.nc mcb750mS2.nc
../refine.py --factor 1 mcb4.nc mcb750mS1.nc
for LEV in 1 2 3 4; do ncdump -h mcb750mS${LEV}.nc |head -n 5; done

# next: chop so have same dims
for LEV in 1 2 3 4; do ncks -O -d x1,,1992 -d y1,,3582 mcb750mS${LEV}.nc mcb750mS${LEV}.nc; done
for LEV in 1 2 3 4; do ncdump -h mcb750mS${LEV}.nc |head -n 5; done

set +x

