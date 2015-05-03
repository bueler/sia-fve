#!/bin/bash

# generate 5 different versions of the 750m MCB bed, with
#   S5 interpolated from the 9000m grid
#   S4 interpolated from the 4500m grid
#   S3 interpolated from the 3000m grid
#   S2 interpolated from the 1500m grid
#   S1 a copy of the 750m grid

set -e # exit on error
set -x

../refine.py --factor 12 mcb0.nc mcb750mS5.nc
../refine.py --factor 6  mcb1.nc mcb750mS4.nc
../refine.py --factor 4  mcb2.nc mcb750mS3.nc
../refine.py --factor 2  mcb3.nc mcb750mS2.nc
../refine.py --factor 1  mcb4.nc mcb750mS1.nc
for LEV in 1 2 3 4 5; do ncdump -h mcb750mS${LEV}.nc |head -n 5; done

# next: chop ocean cells so have same dims
for LEV in 1 2 3 4 5; do ncks -O -d x1,,1980 -d y1,,3576 mcb750mS${LEV}.nc mcb750mS${LEV}.nc; done

# next: copy so all have same cmb and thk; only topg differs
for LEV in 2 3 4 5; do ncks -A -v cmb,thk mcb750mS1.nc mcb750mS${LEV}.nc; done

set +x

