petsc/grn/steadyfinemcb/
========================

After the study in `petsc/grn/robust/`, the purpose of these scripts is to get
a decent high-res solution of the steady-state problem using the MCB bed data.

Run a quickish version
----------------------

First generate the `mcb?.nc` files in `petsc/grn/robust/`.

<!---
# re mcb
#   1=4500m, 2=3000m, 3=1500m, 4=1200m, 5=900m, 6=600m
# re 900 m grid blocks
#   6 | 30, 24, 18, 12, 6
# re 1200 m grid
#   8 | 24, 16, 8
# re 1500 m grid blocks
#   10 | 30, 20, 10
-->
Then do in this directory:

    $ ../refine.py --factor 3 ../robust/mcb1.nc mcb1500mS3.nc
    $ ../refine.py --factor 2 ../robust/mcb2.nc mcb1500mS2.nc
    $ cp ../robust/mcb3.nc mcb1500mS1.nc

Do

    $ ./study.sh  #FIXME

