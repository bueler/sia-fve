petsc/mahaffy/
==============

The code `mahaffy.c` tests Mahaffy and Mahaffy* finite volume element methods,
which solve the shallow ice approximation (SIA) on a structured Q1 grid.  This
code is documented by manuscript `../mahaffy/mahaffyfem.tex`, but also this SIA
calculation is an instance of the general situation documented in
`../paper/lc.tex`.

In any case, `mahaffy.c` uses PETSc DMDA and SNESVI to solve the steady SIA free
boundary problem.

Build the executable, and run on minimal verification case, by

    $ make mahaffy
    $ ./mahaffy

See comments at start of `mahaffy.c` for verification runs.

In addition, `mahaffy.c` can read real data.  In particular, `grn2petsc.py`
reads SeaRISE data from NetCDF format and converts it to PETSc binary format,
for use as input to `mahaffy.c`.  It requires the _netcdf4-python_ module; see
https://github.com/Unidata/netcdf4-python.

Stage 1 of the setup is to get the data as slightly-cleaned NetCDF.  This stage
requires both PISM and NCO.  Note dimensions must be multiples of 3 in
`mahaffy.c` because we want a FD-evaluated Jacobian on a periodic grids with
a stencil width of 1.  Thus we trim the x dimension from 301 values to 300,
but this removes only ocean cells.  The y dimension is already 561, which is
divisible by 3.  Also note 1 / (910 * 31556926) = 3.48228182586954e-11
is the conversion factor to turn  kg m-2 a-1  into  m s-1  for ice with
density 910 kg m-3.

    $ cd ~/pism/examples/std-greenland/
    $ ./preprocess.sh
    $ cd ~/layer-conserve/petsc/mahaffy/    # back to this dir
    $ ln -s ~/pism/examples/std-greenland/pism_Greenland_5km_v1.1.nc
    $ ./cleangrn.sh    # uses NCO to clean up pism_Greenland_5km_v1.1.nc; creates grn.nc

Stage 2 of setup is to convert to PETSc binary using PETSc scripts:

    $ ln -s ~/petsc/bin/petsc-pythonscripts/PetscBinaryIO.py
    $ ln -s ~/petsc/bin/petsc-pythonscripts/petsc_conf.py
    $ python grn2petsc.py grn.nc grn5km.dat    # defaults to 5km

(For maint branch of petsc replace `bin/petsc-pythonscripts/` by just
`bin/petsc-pythonscripts/`).

Stage 3 is to run it and view results:  FIXME: does not make it to level 8; why?

    $ make mahaffy
    $ mkdir test/
    $ ./mahaffy -mah_read grn5km.dat -mah_Neps 8 -mah_showdata -draw_pause 2 -snes_monitor -mah_dump test/
    $ cd test/
    $ ../figsmahaffy.py
    $ eog *.png

A higher-res (2.5 km), parallel, with-upwinding version looks like this; the `-pc_type` and
`-sub_pc_type` choices may be helpful:

    $ python grn2petsc.py --refine 2 grn.nc grn2p5km.dat
    $ mkdir test2p5km/
    $ mpiexec -n 6 ./mahaffy -mah_read grn2p5km.dat -mah_dump test2p5km/ -mah_upwind -snes_monitor -pc_type asm -sub_pc_type lu

Or this

    $ mpiexec -n 6 ./mahaffy -mah_read grn2p5km.dat -mah_upwind -mah_dump test2p5km/ -pc_type asm -sub_pc_type lu -mah_D0 10.0 -snes_max_it 2000
