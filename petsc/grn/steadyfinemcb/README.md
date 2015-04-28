petsc/grn/steadyfinemcb/
========================

After the study in `petsc/grn/robust/`, the purpose of these scripts is to get
a decent high-res solution of the steady-state problem using the MCB bed data.

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

Run the 1500m (quickish) version
--------------------------------

First generate sufficient `mcb?.nc` files in `petsc/grn/robust/` if they are
not already there:

    $ (cd ../robust/ && ./generatemcb.sh 3)

Then do in this directory:

    $ ./setup1500m.sh

To view the resulting files do

    $ ncecat mcb1500mS?.nc foo.nc && ncview foo.nc && rm foo.nc

Do

    $ ./study.sh

