#!/bin/bash

# generate profile comparison figure for bedrock step

set -x

otheropts="-mah_bedstep -mah_D0 0.001 -da_refine 2 -snes_max_it 400 -snes_type vinewtonrsls -pc_type asm -sub_pc_type lu"

./mahaffy -mah_dump $1 $otheropts -mah_noupwind
(cd $1 && mv H.dat Hnoupwind.dat)

./mahaffy -mah_dump $1 $otheropts -mah_upwindfull
(cd $1 && mv H.dat Hupwindfull.dat)

./mahaffy -mah_dump $1 $otheropts

(cd $1 && ../figsmahaffy.py --profile --half --sharpbed -extra_H Hnoupwind.dat,Hupwindfull.dat -extra_H_labels M*,M*noup,M*fullup)

