#!/bin/bash

# do ../quickstart.sh and ./generate.sh first

# FIXME: make optional depending on whether LONG is set

# these correspond to 10km, 5km, 2500m, 1667m, 1000m, 625m grids
#LEVELS="0 1 2 3 5 8"
LEVELS="0 1 2"

set -e # exit on error
set -x

# FIXME: check whether sea${LEV}.nc exist and stop otherwise

if [ -z ${1+x} ]; then
  NN=6
else
  NN=$1
fi

# move old result out of way if present
OUT=study.robust
touch $OUT
mv -f $OUT SAFE_$OUT

# FIXME: loop over sea and mcb

# perform robustness study
STDOPTS="-cs_D0 10.0 -snes_monitor -snes_max_it 200 -pc_type asm -sub_pc_type lu"
for LEV in $LEVELS; do
    DAT=sea${LEV}.dat
    ../../nc2petsc.py sea${LEV}.nc $DAT  # generate .dat; fast
    for METHOD in rsls ssls; do
        DUMP=sea${LEV}_${METHOD}/
        mkdir -p $DUMP
        mpiexec -n $NN ../../mahaffy -mah_read $DAT $STDOPTS -snes_type vinewton${METHOD} -mah_notry -mah_dump $DUMP
        cat $DUMP/history.txt | grep "spacing in x" | sed 's/.* //g' >> $OUT
        cat $DUMP/history.txt | grep "last successful value of eps" | sed 's/.* //g' >> $OUT
        echo $METHOD >> $OUT
    done
done

set +x

echo ""
echo "********  DONE.  SEE RESULT IN $OUT.  ********"

