petsc/grn/robust/
=================

The purpose of these scripts is to study robustness.  The cases studied have:

  * varying bed roughness, because they are based on either SeaRISE (see `../README.md`)
    or mass-conserving-bed (see `../mcb/README.md`) data,

  * varying grid resolution, from 10 km down to 450 m,
  
  * either reduced-set (`-snes_type vinewtonrsls`) or semi-smooth (`vinewtonssls`)
    solution method

FIXME:  Other parameters which could be studied include especially `-cs_D0 XX`.

FIXME:  Study duration `-mah_dtBE` which allows convergence, perhaps only with rsls.

Run the quicker version
-----------------------

First run `quickstart.sh` in `../`.  Then do:

    $ ./generate.sh
    $ ./study.sh       # defaults to NN=6 processes; otherwise run as  ./study.sh NN

Though shortish, this still takes about 20 minutes.  Look at `study.robust`.

Longer version
--------------

On `NN` processors:

    $ LONG=1 ./generate.sh NN
    $ LONG=1 ./study.sh NN

Generate figure in paper
------------------------

FIXME

Clean up the mess
-----------------

    $ make clean

