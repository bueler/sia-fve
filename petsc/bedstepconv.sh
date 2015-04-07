#!/bin/bash

# refinement path for bedrock step test
# for -da_refine levels 0 ... 3 in file foo.verif, on 6 processes, do
#   $ ./bedstepconv.sh foo.verif 3 6

# creates and destroys files lev?history.txt

for lev in `seq 0 $2`; do
   rm -f lev${lev}history.txt
   cmd="mpiexec -n $3 ./mahaffy -mah_bedstep -mah_D0 0.001 -da_refine $lev -snes_max_it 400 -mah_history lev${lev} -pc_type asm -sub_pc_type lu -snes_type vinewtonrsls -mah_Neps 12 $4"
   echo $cmd
   $cmd
done

rm -f $1
echo $(($lev+1)) > $1
cat lev?history.txt | grep "spacing in x"        | sed 's/.* //g' >> $1
cat lev?history.txt | grep "last successful value of eps"     | sed 's/.* //g' >> $1
cat lev?history.txt | grep "max thickness error" | sed 's/.* //g' >> $1
cat lev?history.txt | grep "av thickness error"  | sed 's/.* //g' >> $1
cat lev?history.txt | grep "solution ice volume" | sed 's/.* //g' >> $1
cat lev?history.txt | grep "exact ice volume"    | sed 's/.* //g' >> $1

rm -f lev?history.txt

