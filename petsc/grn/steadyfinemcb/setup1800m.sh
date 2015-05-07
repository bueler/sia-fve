#!/bin/bash

# generate 4 different versions of the 1800m bed, with
#   S4 interpolated from the 9000m grid mbc0.nc
#   S3 interpolated from the 5400m grid mbc1.nc
#   S2 interpolated from the 3600m grid mbc2.nc
#   S1 a copy of         the 1800m grid mcb3.nc

set -e # exit on error
set -x

../refine.py --factor 5 mcb0.nc mcb1800mS4.nc
../refine.py --factor 3 mcb1.nc mcb1800mS3.nc
../refine.py --factor 2 mcb2.nc mcb1800mS2.nc
../refine.py --factor 1 mcb3.nc mcb1800mS1.nc
for LEV in 1 2 3 4; do ncdump -h mcb1800mS${LEV}.nc |head -n 5; done

# next: chop so have same dims
for LEV in 1 2 3 4; do ncks -O -d x1,,825 -d y1,,1490 mcb1800mS${LEV}.nc mcb1800mS${LEV}.nc; done
for LEV in 1 2 3 4; do ncdump -h mcb1800mS${LEV}.nc |head -n 5; done

set +x

