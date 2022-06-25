#!/bin/bash

# measure serial performance scaling refinement path for dome test
# for -da_refine levels 0 ... 4 in file foo.perf, on 1 processes, do
#   $ ./domeperf.sh foo.perf 4

# creates and destroys files lev?history.txt

rm -f $1
for lev in `seq 0 $2`; do
   cmd="./mahaffy -da_refine $lev -cs_D0 1.0 -snes_type vinewtonssls -snes_max_it 400 -pc_type lu -mah_notry -log_view"
   echo $cmd
   rm -f lev${lev}history.txt
   $cmd &> lev${lev}history.txt
   grep "grid of " lev${lev}history.txt >> $1
   grep "Flop:      " lev${lev}history.txt >> $1
done
