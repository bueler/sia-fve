petsc/grn/phoenix/
==================

The purpose of these scripts is to get a decent high-res solution of the
steady-state problem using the MCB bed data.

The approaches in `petsc/grn/robust/` and `petsc/grn/steadyfinemcb/` both have
problems.
 
Coarsen the 150m MCB data
-------------------------

First generate sufficient `mcb?.nc` files if they are not already present (takes
about 30 minutes):

    $ ln -s ../mcb/MCdataset-2014-11-19.nc    # or download with ../mcb/getmcb.sh
    $ ./generate.sh

Run
---

Do this:

    $ ./study.sh NN

