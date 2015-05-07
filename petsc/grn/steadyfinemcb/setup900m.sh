#!/bin/bash

# generate 5 different versions of the 900m bed, with
#   S5 interpolated from the 9000m grid mbc0.nc
#   S4 interpolated from the 5400m grid mbc1.nc
#   S3 interpolated from the 3600m grid mbc2.nc
#   S2 interpolated from the 1800m grid mbc3.nc
#   S1 a copy of         the  900m grid mcb4.nc

set -e # exit on error
set -x

../refine.py --factor 10 mcb0.nc mcb900mS5.nc
../refine.py --factor 6  mcb1.nc mcb900mS4.nc
../refine.py --factor 4  mcb2.nc mcb900mS3.nc
../refine.py --factor 2  mcb3.nc mcb900mS2.nc
../refine.py --factor 1  mcb4.nc mcb900mS1.nc
for LEV in 1 2 3 4 5; do ncdump -h mcb900mS${LEV}.nc |head -n 5; done

# next: chop somewhat symmetrically so have same dims
ncks -O -d x1,,1650  -d y1,,2980  mcb900mS5.nc mcb900mS5.nc  # 1650-0=1650, 2980-0=2980
ncks -O -d x1,6,1656 -d y1,1,2981 mcb900mS4.nc mcb900mS4.nc  # 1656-6=1650, 2981-1=2980
ncks -O -d x1,7,1657 -d y1,2,2982 mcb900mS3.nc mcb900mS3.nc  # 1657-7=1650, 2982-2=2980
ncks -O -d x1,8,1658 -d y1,4,2984 mcb900mS2.nc mcb900mS2.nc  # 1658-8=1650, 2984-4=2980
ncks -O -d x1,9,1659 -d y1,5,2985 mcb900mS1.nc mcb900mS1.nc  # 1659-9=1650, 2985-5=2980

for LEV in 1 2 3 4 5; do ncdump -h mcb900mS${LEV}.nc |head -n 5; done

set +x

