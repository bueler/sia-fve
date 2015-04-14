#!/bin/bash

set -e # exit on error
set -x # echo commands

ncks -O -v x1,y1,thk,topg,climatic_mass_balance pism_Greenland_5km_v1.1.nc grn.nc

./inplace.py --fixcmbunits --ranges grn.nc

./inplace.py --oceanfix --ranges grn.nc

#./inplace.py --bedundip --sweeps 2 --ranges grn.nc

./nc2petsc.py grn.nc grn.dat

(cd ../ && make mahaffy)

mkdir -p test/

mpiexec -n 6 ../mahaffy -mah_read grn.dat -mah_showdata -draw_pause 2 -snes_monitor -mah_dump test/

