#!/bin/bash

# refinement path for dome test: for result from -da_refine levels 0 ... 4
# put into file foo.verif, on 6 processes, do
#   $ ./domeconv.sh foo.verif 4 6

# creates and destroys files lev?history.txt

for lev in `seq 0 $2`; do
   rm -f lev${lev}history.txt
   mpiexec -n $3 ./mahaffy -da_refine $lev -snes_max_it 200 -pc_type asm -sub_pc_type lu -mah_history lev${lev}
done

rm -f $1
echo $(($lev+1)) > $1
cat lev?history.txt | grep "spacing in x"        | sed 's/.* //g' >> $1
cat lev?history.txt | grep "last successful"     | sed 's/.* //g' >> $1
cat lev?history.txt | grep "max thickness error" | sed 's/.* //g' >> $1
cat lev?history.txt | grep "av thickness error"  | sed 's/.* //g' >> $1

rm -f lev?history.txt

