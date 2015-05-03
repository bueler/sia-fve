petsc/grn/steadyfinemcb/
========================

After the study in `petsc/grn/robust/`, the purpose of these scripts is to get
a decent high-res solution of the steady-state problem using the MCB bed data.

Coarsen the 150m MCB data
-------------------------

<!---
# re mcb
#   1=4500m, 2=3000m, 3=1500m, 4=750m
-->
First generate sufficient `mcb?.nc` files if they are not already present (takes
about 30 minutes):

    $ ln -s ../mcb/MCdataset-2014-11-19.nc    # or download with ../mcb/getmcb.sh
    $ ./generate.sh

Run the quicker 1500m version
-----------------------------

<!---
# re 1500 m grid blocks
#   10 | 30, 20, 10
-->
Do this:

    $ ./setup1500m.sh

To view the resulting files do something like this:

    $ ncecat mcb1500mS?.nc foo.nc && ncview foo.nc && rm foo.nc

This version requires about 4 Gb memory.  Do the runs on the default NN=6
MPI processes:

    $ ./study1500m.sh

Do this for NN processes:

    $ ./study1500m.sh NN


Run the 750m version
-----------------------------

_Warning: requires at least 12 Gb and probably more._

<!---
# re 750 m grid blocks
#   5 | 60, 30, 20, 10, 5
-->

Setup and view:

    $ ./setup750m.sh    # maybe 60 minutes?
    $ ncecat mcb750mS?.nc foo.nc && ncview foo.nc && rm foo.nc

Do this for NN processors:

    $ ./study750m.sh NN
