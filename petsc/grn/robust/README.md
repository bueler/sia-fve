petsc/grn/robust/
=================

The purpose of these scripts is to study robustness.  The cases studied have:

  * varying bed roughness, because they are based on either SeaRISE (see `../README.md`)
    or mass-conserving-bed (=MCB; see `../mcb/README.md`) data,

  * varying grid resolution, from 10 km down to 450 m,
  
  * either reduced-set (`-snes_type vinewtonrsls`) or semi-smooth (`vinewtonssls`)
    solution method

Run a quickish version
----------------------

First run `quickstart.sh` in the `grn/` directory and the download script in the
`grn/mcb/` directory:

    $ (cd ../ && QSSHOW= ./quickstart.sh)
    $ (cd ../mcb/ && ./getmcb.sh)     # downloads with wget; no action if present
    $ ln -s ../mcb/MCdataset-2014-11-19.nc

Then do:

    $ ./generatesea.sh 2

to generate files `sea?.nc` and

    $ ./generatemcb.sh 2

to generate files `mcb?.nc`.  In the latter you may want to set the PISM
location, e.g.

    $ PISM=/u1/uaf/bueler/pism ./generatemcb.sh 2

for me on pacman.

The above generates down to 2500m grid for SeaRISE and 4500m for MCB.  Then do:

    $ ./study.sh

The above defaults to NN=6 processes.  For different number NN of processes, do

    $ ./study.sh NN

Though shortish, this still takes about 20 minutes.  Look at `study.robust`.

Longest version
---------------

This goes all the way to 600m grids on both SeaRISE and MCB sequences.
On `NN` processors:

    $ ./generatesea.sh 6
    $ PISM=/u1/uaf/bueler/pism ./generatemcb.sh 6
    $ ./study.sh NN

Generate figure in paper
------------------------

Copy the result `study.robust` to `../../../paper/`.

Clean up the mess
-----------------

    $ make clean

