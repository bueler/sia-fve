#!/bin/bash

for NN in 1 2 4 ; do
  for sch in 0 1 2 ; do
    echo "testing advection scheme $sch on $NN processes:" ;
    for lev in 0 2 4 6 8 10 12 ; do
      mpiexec -n $NN ./layer -lay_exactinit -lay_noshow -lay_dt 0.01 -lay_steps 10 -lay_adscheme $sch -da_refine $lev -lay_jac | 'grep' error
    done
    echo
  done
done
