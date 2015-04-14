#!/bin/bash

set -e # exit on error
set -x # echo commands

ncks -O -v x1,y1,thk,topg,climatic_mass_balance pism_Greenland_5km_v1.1.nc grn.nc

./inplace.py --fixcmbunits --ranges grn.nc

./inplace.py --oceanfix --ranges grn.nc

### OPTIONAL ###

# refine:
#./refine.py --factor 2 grn.nc grn2p5km.nc

# coarsen:
#ncks -d x1,,,4 -d y1,,,4 grn.nc grn20km.nc

# apply bed sweeps:
#cp grn.nc grn_origtopg.nc
#./inplace.py --bedundip --sweeps 2 --ranges grn.nc

# trim so that -snes_fd_color works:
#ncks -O -d x1,0,299 grn.nc grn.nc

./nc2petsc.py grn.nc grn.dat

mkdir -p test/

mpiexec -n 6 ../mahaffy -mah_read grn.dat -mah_D0 1.0 -mah_showdata -draw_pause 2 -snes_monitor -mah_dump test/ -pc_type asm -sub_pc_type lu -mah_notry

### RUN OPTIONS ###
#--mah_dtBE 10.0  <-- replaces -mah_notry

