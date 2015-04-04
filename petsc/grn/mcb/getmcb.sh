#!/bin/bash

wget -nc ftp://sidads.colorado.edu/DATASETS/IDBMG4_BedMachineGr/MCdataset-2014-11-19.nc

#ncks -v x,y,thickness,bed MCdataset-2014-11-19.nc -d x,,,6 -d y,,,6 mcb900m.nc
ncks -v x,y,thickness,bed MCdataset-2014-11-19.nc -d x,,,30 -d y,,,30 mcb4500m.nc

