#!/bin/bash

DATANAME=MCdataset-2014-11-19.nc
echo "downloading 1.7Gb file $DATANAME with 150 m grid ..."
wget -nc ftp://sidads.colorado.edu/DATASETS/IDBMG4_BedMachineGr/$DATANAME

