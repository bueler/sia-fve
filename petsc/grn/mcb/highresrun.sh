#!/bin/bash

# use like this:
#   ./highresrun.sh 30 4500m 6
# to get 4500m run on 6 processes

set -e # exit on error
set -x # echo commands

mkdir test$2/     # immediately generates error if already exists

./getmcb.sh
./average.py --block $1 mcb$2.nc
./remap2mcb.py searise5km.nc mcb$2.nc fix$2.nc
../grn2petsc.py fix$2.nc fix$2.dat

mpiexec -n $3 ../../mahaffy -mah_read fix$2.dat -snes_monitor -pc_type asm -sub_pc_type lu -snes_max_it 200 -mah_dump test$2/ -mah_dtBE 10.0
