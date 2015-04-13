#!/bin/bash

ncks -O -v x1,y1,thk,topg,climatic_mass_balance pism_Greenland_5km_v1.1.nc -d x1,0,299 grn.nc
# convert  kg m-2 year-1  to  m s-1  for ice of 910.0 kg m-3
ncap2 -O -s "cmb=3.48228182586954e-11*climatic_mass_balance" grn.nc grn.nc
ncatted -O -a units,cmb,o,c,"m s-1" grn.nc

