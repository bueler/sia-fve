#!/bin/bash

SKIP=$1
OUTNAME=$2

DATANAME=MCdataset-2014-11-19.nc
echo "downloading $DATANAME with 150 m grid ..."
wget -nc ftp://sidads.colorado.edu/DATASETS/IDBMG4_BedMachineGr/$DATANAME

echo "generating $OUTNAME by stride of $SKIP ..."
ncks -v x,y,thickness,bed MCdataset-2014-11-19.nc -d x,,,$SKIP -d y,,,$SKIP $OUTNAME
ncrename -O -v x,x1 $OUTNAME
ncrename -O -v y,y1 $OUTNAME
ncrename -O -v thickness,thk $OUTNAME
ncrename -O -v bed,topg_nobathy $OUTNAME
