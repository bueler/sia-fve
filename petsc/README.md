layer-conserve/petsc/
==============

layer.c
-------

The code `layer.c` describes a one-dimensional moving layer with a non-negative-
constrained thickness, and thus a moving boundary.  PETSc SNESVI is used to
solve the free boundary problem at each time step.  For more information build
`doclayer.pdf` by

    $ make doclayer.pdf

Building `layer.c` requires PETSc 3.5.2 or maint.  Do

    $ make layer

To thoroughly test with exact solution,

    $ ./convtest.sh


mahaffy.c
---------

Tests Mahaffy and Mahaffy* finite volume element methods, which solve the
shallow ice approximation (SIA) on a structured Q1 grid.  Documented by manuscript
`../mahaffy/mahaffyfem.tex`.  This SIA calculation is an instance of the
general situation documented in `../paper/lc.tex`, and again PETSc SNESVI is used to
solve the free boundary problem.

Build the executable, and run on minimal verification case, by

    $ make mahaffy
    $ ./mahaffy

See comments at start of `mahaffy.c` for runs.


grn2petsc.py
------------

Read SeaRISE data from NetCDF format and convert to PETSc binary format, for
use as input to `mahaffy.c`.  Requires the _netcdf4-python_ module; see
https://github.com/Unidata/netcdf4-python.   some setup:

    $ cd ~/pism/examples/std-greenland/
    $ ./preprocess.sh
    $ cd ~/layer-conserve/petsc/  # back to this dir
    $ ln -s ~/pism/examples/std-greenland/pism_Greenland_5km_v1.1.nc
    $ ln -s ~/petsc-maint/bin/pythonscripts/PetscBinaryIO.py
    $ ln -s ~/petsc-maint/bin/pythonscripts/petsc_conf.py
    $ python grn2petsc.py pism_Greenland_5km_v1.1.nc grn.dat
    $ make mahaffy
    $ ./mahaffy -mah_i grn.dat FIXME

