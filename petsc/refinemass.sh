#!/bin/bash

NN=1

for MM in 0 1 2 ; do
  STEPS=$(bc -l <<< "100*10^$MM")
  DT=$(bc -l <<< "scale = 6; 10.0/$STEPS")
  cmd="mpiexec -n $NN ./layer -da_refine 5 -lay_dt ${DT} -lay_steps ${STEPS} -lay_adscheme 1 -lay_jac -lay_noshow -lay_timedependentsource -lay_massfile mass${MM}.txt -snes_max_it 1000"
  echo $cmd
  $cmd
done

# now do something like
#   ./masstimefig.py -o foo.pdf mass?.txt

