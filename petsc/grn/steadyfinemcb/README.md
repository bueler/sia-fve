petsc/grn/steadyfinemcb/
========================

After the study in `petsc/grn/robust/`, the purpose of these scripts is to get
a decent high-res solution of the steady-state problem using the MCB bed data.

First generate sufficient `mcb?.nc` files if they are not already present (takes
about 30 minutes):

    $ ln -s ../mcb/MCdataset-2014-11-19.nc    # or download with ../mcb/getmcb.sh
    $ ./generate.sh


Run the quicker 1800m version
-----------------------------

Do this:

    $ ./setup1800m.sh    # maybe 15 minutes?

To view the resulting files do something like this:

    $ ncecat mcb1800mS?.nc foo.nc && ncview foo.nc && rm foo.nc

This version requires about 4 Gb memory.  Do the runs on the default NN=6
MPI processes:

    $ ./study.sh

Do this for NN processes:

    $ ./study.sh NN


Run the 900m version
-----------------------------

Setup and view:

    $ ./setup900m.sh    # maybe 60 minutes?
    $ ncecat mcb900mS?.nc foo.nc && ncview foo.nc && rm foo.nc

Do this for NN processors:

    $ FINE=1 ./study.sh NN

