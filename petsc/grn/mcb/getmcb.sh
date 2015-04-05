#!/bin/bash

DATANAME=MCdataset-2014-11-19.nc
echo "downloading $DATANAME with 150 m grid ..."
wget -nc ftp://sidads.colorado.edu/DATASETS/IDBMG4_BedMachineGr/$DATANAME

