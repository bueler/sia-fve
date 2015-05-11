#!/bin/bash

# an example short run on 4 cores for testing:
#   RUNT=0.01 RUNDT=0.001 RUNDUMPDT=0.0025 ./keepgoing.sh 4 unnamed.dat bar/

set -e # exit on error
set -x

NN=$1

if [ -z ${2+x} ]; then
  echo "only runs with two or three arguments:  keepgoing.sh NN FILE [DUMPDIR/]"
  exit
else
  INFILE=$2
fi

if [ -z ${3+x} ]; then
  DUMP=dump/
else
  DUMP=$3
fi

if [ -z ${RUNT+x} ]; then
  RUNT=1.0
fi

if [ -z ${RUNDT+x} ]; then
  RUNDT=0.05
fi

if [ -z ${RUNDUMPDT+x} ]; then
  RUNDUMPDT=-1.0  # inactive when this is negative
fi

EXEC=../../mahaffy
STDOPTS="-cs_D0 30.0 -snes_monitor -pc_type asm -sub_pc_type lu"

mkdir -p $DUMP
mpiexec -n $NN $EXEC -mah_read $INFILE -mah_readinitialsurface $INFILE \
        $STDOPTS -mah_T $RUNT -mah_dt $RUNDT -cs_start 12 \
        -mah_dumpdt $RUNDUMPDT -mah_dump $DUMP

set +x

echo ""
echo "********  DONE.  ********"

