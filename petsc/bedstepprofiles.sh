#!/bin/bash

# generate profile comparison figure for bedrock step

set -x

otheropts="-mah_bedstep -mah_D0 0.001 -da_refine 1 -snes_max_it 400 -pc_type asm -sub_pc_type lu -mah_Neps 12"

./mahaffy -mah_dump $1 $otheropts -mah_lambda 0.0
(cd $1 && mv H.dat Hnoupwind.dat)

./mahaffy -mah_dump $1 $otheropts -mah_lambda 1.0
(cd $1 && mv H.dat Hupwindfull.dat)

./mahaffy -mah_dump $1 $otheropts

(cd $1 && ../figsmahaffy.py --profile --half --exactbed -extra_H Hnoupwind.dat,Hupwindfull.dat -extra_H_labels M*,M*no,M*full)

