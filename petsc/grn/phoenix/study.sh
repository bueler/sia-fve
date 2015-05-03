#!/bin/bash

set -e # exit on error
set -x

if [ -z ${1+x} ]; then
  NN=6
else
  NN=$1
fi

# move old result out of way if present
OUT=study.750m
touch $OUT
mv -f $OUT SAFE_$OUT

STDOPTS="-cs_D0 10.0 -snes_monitor -snes_max_it 400 -pc_type asm -sub_pc_type lu"

extracthistory ()
{
cat $1/history.txt | grep "spacing in x-direction" | sed 's/.* //g' >> $2
cat $1/history.txt | grep "last successful value of eps" | sed 's/.* //g' >> $2
cat $1/history.txt | grep "maximum solution diffusivity" | sed 's/.* //g' >> $2
cat $1/history.txt | grep "total time" | sed 's/.* //g' >> $2
}

NAMEROOT=mcb
DAT=${NAMEROOT}4.dat
../../nc2petsc.py ${NAMEROOT}4.nc $DAT
DUMP=m4/
mkdir -p $DUMP
mpiexec -n $NN ../../mahaffy -mah_read $DAT $STDOPTS -mah_notry -mah_dump $DUMP
head -n 1 $DUMP/history.txt >> $OUT
extracthistory $DUMP $OUT

LEVELS="3 2 1 0"
PREV=$DUMP/unnamed.dat
for LEV in $LEVELS; do
    DAT=${NAMEROOT}${LEV}.dat
    # refine previous result
    ../../petsc2nc.py $PREV tmp.nc
    FIXME ../refine.py --factor 2 tmp.nc ${NAMEROOT}${LEV}.nc
    # copy fine bed into refined-previous
    FIXME ncks -A -v cmb,thk A.nc ${NAMEROOT}${LEV}.nc
    ../../nc2petsc.py ${NAMEROOT}${LEV}.nc $DAT
    DUMP=m${LEV}/
    mkdir -p $DUMP
    mpiexec -n $NN ../../mahaffy -cs_start 8 -mah_read $DAT -mah_readinitial $DAT $STDOPTS -mah_dump $DUMP
    head -n 1 $DUMP/history.txt >> $OUT
    extracthistory $DUMP $OUT
done

set +x

echo ""
echo "********  DONE.  SEE RESULT IN $OUT.  ********"

