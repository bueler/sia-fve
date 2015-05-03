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
cat ${1}history.txt | grep "spacing in x-direction" | sed 's/.* //g' >> $2
cat ${1}history.txt | grep "last successful value of eps" | sed 's/.* //g' >> $2
cat ${1}history.txt | grep "maximum solution diffusivity" | sed 's/.* //g' >> $2
cat ${1}history.txt | grep "total time" | sed 's/.* //g' >> $2
}

NAMEROOT=mcb

NC=${NAMEROOT}1.nc
DAT=${NAMEROOT}1.dat
DUMP=m1/

# merge refined topg from coarser level
../refine.py --factor 2 ${NAMEROOT}0.nc tmptopgrefine.nc
cp $NC tmp${NC}
ncks -A -v topg tmptopgrefine.nc tmp${NC}

# generate .dat for run
../../nc2petsc.py tmp${NC} $DAT
rm tmptopgrefine.nc tmp${NC}
mkdir -p $DUMP

# run twice
mpiexec -n $NN ../../mahaffy -mah_read $DAT $STDOPTS -mah_dtrecovery 10.0 -mah_dump $DUMP
cp ${DUMP}unnamed.dat foo.dat
rm -rf $DUMP
mkdir -p $DUMP
mpiexec -n $NN ../../mahaffy -mah_read foo.dat -mah_readinitial foo.dat -cs_start 9 $STDOPTS -mah_dtrecovery 10.0 -mah_dump $DUMP
cp ${DUMP}unnamed.dat foo.dat
rm -rf $DUMP
mkdir -p $DUMP
mpiexec -n $NN ../../mahaffy -mah_read foo.dat -mah_readinitial foo.dat -cs_start 10 $STDOPTS -mah_dtrecovery 10.0 -mah_dump $DUMP

head -n 1 ${DUMP}history.txt >> $OUT
extracthistory $DUMP $OUT

LEVELS="2"
#LEVELS="2 3 4"
PREVLEV=1
for LEV in $LEVELS; do
    PREV=${DUMP}unnamed.dat
    NC=${NAMEROOT}${LEV}.nc
    DAT=${NAMEROOT}${LEV}.dat
    DUMP=m${LEV}/

    # merge refined topg from coarser level
    ../refine.py --factor 2 ${NAMEROOT}${PREVLEV}.nc tmptopgrefine.nc
    cp $NC tmp${NC}
    ncks -A -v topg tmptopgrefine.nc tmp${NC}

    # merge refined previous result for thk
    ../../petsc2nc.py $PREV tmp.nc
    ../refine.py --factor 2 tmp.nc tmpthkrefine.nc
    ncks -A -v thk tmpthkrefine.nc tmp${NC}

    # generate .dat for run
    ../../nc2petsc.py tmp${NC} $DAT
    rm tmptopgrefine.nc tmpthkrefine.nc tmp.nc tmp${NC}
    mkdir -p $DUMP

    # FIXME: next line *could* be time-stepping, but it needs to do a lot of it
    mpiexec -n $NN ../../mahaffy -cs_start 5 -mah_read $DAT -mah_readinitial $DAT $STDOPTS -mah_dtrecovery 10.0 -mah_dump $DUMP
    mpiexec -n $NN ../../mahaffy -cs_start 7 -mah_read ${DUMP}unnamed.dat -mah_readinitial ${DUMP}unnamed.dat $STDOPTS -mah_dtrecovery 10.0 -mah_dump $DUMP

    head -n 1 ${DUMP}history.txt >> $OUT
    extracthistory $DUMP $OUT

    PREVLEV=$LEV
done

set +x

echo ""
echo "********  DONE.  SEE RESULT IN $OUT.  ********"

