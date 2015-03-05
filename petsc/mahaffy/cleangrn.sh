#!/bin/bash

ncks -v x1,y1,thk,topg,climatic_mass_balance pism_Greenland_5km_v1.1.nc -d x1,0,299 grn.nc
ncap2 -O -s "cmb=3.48228182586954e-11*climatic_mass_balance" grn.nc grn.nc
ncatted -O -a units,cmb,o,c,"m s-1" grn.nc

