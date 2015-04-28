#!/bin/bash

# generate 3 different versions of the 1500m MCB bed, with S3 smoothed from the
# 4500m grid and S2 smoothed from the 3000m grid

set -e # exit on error
set -x

../refine.py --factor 3 ../robust/mcb1.nc mcb1500mS3.nc
../refine.py --factor 2 ../robust/mcb2.nc mcb1500mS2.nc
../refine.py --factor 1 ../robust/mcb3.nc mcb1500mS1.nc
for LEV in 1 2 3; do ncdump -h mcb1500mS${LEV}.nc |head -n 5; done    # next: chop so have same dims
for LEV in 1 2 3; do ncks -O -d x1,,996 -d y1,,1791 mcb1500mS${LEV}.nc mcb1500mS${LEV}.nc; done
for LEV in 1 2 3; do ncdump -h mcb1500mS${LEV}.nc |head -n 5; done

set +x

