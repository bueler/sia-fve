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
https://github.com/Unidata/netcdf4-python.

Stage 1 of setup is to get the data as slightly-cleaned NetCDF.  This stage
requires both PISM and NCO.  Note dimensions must be multiples of 3 in
`mahaffy.c` because we want a FD-evaluated Jacobian on a periodic grids with
a stencil width of 1.  Thus we trim the x dimension from 301 values to 300,
but this removes only ocean cells.  The y dimension is already 561, which is
divisible by 3.  Also note 1 / (910 * 31556926) = 3.48228182586954e-11
is the conversion factor to turn  kg m-2 a-1  into  m s-1  for ice with
density 910 kg m-3.

    $ cd ~/pism/examples/std-greenland/
    $ ./preprocess.sh
    $ cd ~/layer-conserve/petsc/  # back to this dir
    $ ln -s ~/pism/examples/std-greenland/pism_Greenland_5km_v1.1.nc
    $ ncks -v x1,y1,thk,topg,climatic_mass_balance pism_Greenland_5km_v1.1.nc -d x1,0,299 grn.nc
    $ ncap2 -O -s "cmb=3.48228182586954e-11*climatic_mass_balance" grn.nc grn.nc
    $ ncatted -O -a units,cmb,o,c,"m s-1" grn.nc

Stage 2 of setup is to convert to PETSc binary using PETSc scripts:

    $ ln -s ~/petsc-maint/bin/pythonscripts/PetscBinaryIO.py
    $ ln -s ~/petsc-maint/bin/pythonscripts/petsc_conf.py
    $ python grn2petsc.py grn.nc grn.dat

Stage 3 is to run it and view results:

    $ make mahaffy
    $ mkdir test/
    $ ./mahaffy -mah_read grn.dat -mah_Neps 8 -mah_showdata -draw_pause 2 -snes_monitor -mah_dump test/
    $ cd test/
    $ ../figsmahaffy.py
    $ eog *.png

FIXME: above run does not make it to level 8

FIXME: seems not to run in parallel:
    $ mpiexec -n 6 ./mahaffy -mah_read grn.dat -mah_Neps 8 -snes_max_it 200 -snes_monitor -pc_type asm -sub_pc_type lu

