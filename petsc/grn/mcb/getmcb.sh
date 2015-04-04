#!/bin/bash

DATANAME=MCdataset-2014-11-19.nc
echo "downloading $DATANAME ..."
wget -nc ftp://sidads.colorado.edu/DATASETS/IDBMG4_BedMachineGr/$DATANAME

COARSE=mcb4500m.nc
echo "generating $COARSE ..."
#ncks -v x,y,thickness,bed MCdataset-2014-11-19.nc -d x,,,6 -d y,,,6 mcb900m.nc
ncks -v x,y,thickness,bed MCdataset-2014-11-19.nc -d x,,,30 -d y,,,30 $COARSE
ncrename -O -v x,x1 $COARSE
ncrename -O -v y,y1 $COARSE
ncrename -O -v thickness,thk $COARSE
ncrename -O -v bed,topg_nobathy $COARSE
