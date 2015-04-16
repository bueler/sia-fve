#!/bin/bash

# refinement path for dome test
# for -da_refine levels 0 ... 4 in file foo.verif, on 6 processes, do
#   $ ./domeconv.sh foo.verif 4 6
# for same plus -mah_true:
#   $ ./domeconv.sh foo.verif 4 6 -mah_true

# creates and destroys files lev?history.txt

for lev in `seq 0 $2`; do
   rm -f lev${lev}history.txt
   cmd="mpiexec -n $3 ./mahaffy $4 -da_refine $lev -cs_D0 1.0 -snes_type vinewtonssls -snes_max_it 400 -pc_type asm -sub_pc_type lu -mah_history lev${lev} -mah_notry"
   echo $cmd
   $cmd
done

rm -f $1
echo $(($lev+1)) > $1
cat lev?history.txt | grep "spacing in x"        | sed 's/.* //g' >> $1
cat lev?history.txt | grep "last successful value of eps"     | sed 's/.* //g' >> $1
cat lev?history.txt | grep "max thickness error" | sed 's/.* //g' >> $1
cat lev?history.txt | grep "av thickness error"  | sed 's/.* //g' >> $1

rm -f lev?history.txt

