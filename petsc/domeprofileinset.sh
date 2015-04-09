#!/bin/bash

# generate profile comparison figure for bedrock step

set -x

otheropts="-mah_D0 1.0 -snes_type vinewtonssls -snes_max_it 100"

./mahaffy -mah_dump $1 $otheropts -da_refine 3
(cd $1 && mv x.dat xfiner.dat)
(cd $1 && mv y.dat yfiner.dat)
(cd $1 && mv H.dat Hfiner.dat)

./mahaffy -mah_dump $1 $otheropts -da_refine 2

(cd $1 && ../figsmahaffy.py --profile --half --exactdome --blowup -extra_H Hfiner.dat -extra_H_labels M*,foo)

