#!/bin/bash

# do ./generatesea.sh and ./generatemcb.sh first

set -e # exit on error
set -x

if [ -z ${1+x} ]; then
  NN=6
else
  NN=$1
fi

# move old result out of way if present
OUT=study.robust
touch $OUT
mv -f $OUT SAFE_$OUT

# perform robustness study
LEVELS="0 1 2 3 4 5 6"
STDOPTS="-cs_D0 10.0 -snes_monitor -snes_max_it 400 -pc_type asm -sub_pc_type lu"
for GRN in sea mcb; do
    for LEV in $LEVELS; do
        STEM=${GRN}${LEV}
        NC=${STEM}.nc
        if [ -e $NC ]; then
            DAT=${STEM}.dat
            ../../nc2petsc.py $NC $DAT  # generate .dat; fast
            for METHOD in rsls ssls; do
                DUMP=${STEM}_${METHOD}/
                mkdir -p $DUMP
                mpiexec -n $NN ../../mahaffy -mah_read $DAT $STDOPTS -snes_type vinewton${METHOD} -mah_notry -mah_dump $DUMP
                cat $DUMP/history.txt | grep "spacing in x" | sed 's/.* //g' >> $OUT
                cat $DUMP/history.txt | grep "last successful value of eps" | sed 's/.* //g' >> $OUT
                cat $DUMP/history.txt | grep "total time" | sed 's/.* //g' >> $OUT
                echo $GRN >> $OUT
                echo $METHOD >> $OUT
            done
        fi
    done
done

set +x

echo ""
echo "********  DONE.  SEE RESULT IN $OUT.  ********"

