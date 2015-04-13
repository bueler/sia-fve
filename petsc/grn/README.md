petsc/grn/
==========

The main script `grn2petsc.py` here reads SeaRISE data from NetCDF format and
converts it to PETSc binary format.  It requires the
[netcdf4-python module](https://github.com/Unidata/netcdf4-python).
The code `mahaffy` can read the specially-formatted binary file it produces.

Get the data from the PISM source
---------------------------------

This stage requires [NCO](http://nco.sourceforge.net/) and a
copy of [PISM](http://www.pism-docs.org).  Do:

    $ cd ~/pism/examples/std-greenland/
    $ ./preprocess.sh                   # download Greenland_5km_v1.1.nc with wget if not present
    $ cd ~/sia-fve/petsc/grn/           # back to this dir
    $ ln -s ~/pism/examples/std-greenland/pism_Greenland_5km_v1.1.nc

Clean the data
--------------

These stages all work on the data as NetCDF files, so their results can be
viewed with `ncview` or similar.

The dimensions in the `.nc` file must be multiples of 3 if we want `mahaffy.c`
to use a FD-by-coloring-evaluated Jacobian on a periodic
grid with stencil width 1.  Thus we trim the x dimension from 301 values to 300,
but this removes only ocean cells.  The y dimension is already 561, which is
divisible by 3.

    $ ./cleangrn.sh         # creates grn.nc by trim dimensions and removing unwanted vars

Also we convert the climatic mass balance variable kg m-2 a-1  into  m s-1:

    $ ./inplace.py --fixcmbunits --ranges grn.nc

Also we change the ocean topography to not be so deep, because the ice flow
model does not know about flotation.  And we change the ocean value of climatic
mass balance to very negative, because this simulates strong calving.

    $ ./inplace.py --oceanfix --ranges grn.nc

FIXME: not implemented:
Optionally we apply sweeps to remove bed dips:

    $ ./inplace.py --bedundip --sweeps 2 --ranges grn.nc

Convert to PETSC binary
-----------------------

This stage uses PETSc scripts, with optional refinement.

    $ ln -s ~/petsc/bin/petsc-pythonscripts/PetscBinaryIO.py  # or similar for maint/3.5
    $ ln -s ~/petsc/bin/petsc-pythonscripts/petsc_conf.py
    $ ./nc2petsc.py grn.nc grn5km.dat                         # defaults to 5km

If desired one can go to a finer grid by bilinear interpolation; this uses a
module called `interpad.py`:

    $ ./nc2petsc.py --refine 5 grn.nc grn1km.dat

Run the model
-------------

Build the executable `mahaffy`, and run on 5km grid using an over-estimate of
diffusivity:

    $ (cd ../ && make mahaffy)
    $ mkdir test/
    $ mpiexec -n 6 ../mahaffy -mah_read grn5km.dat -mah_showdata -draw_pause 2 -snes_monitor -mah_dump test/

This run only takes a few minutes and uses the data as is.

Generate figures
----------------

Generate `.pdf` and `.png` figures:

    $ cd test/
    $ ../../figsmahaffy.py --observed

Higher/lower resolution
-----------------------

A higher-res (2.5 km), parallel version might look like this; the robust
`-pc_type asm -sub_pc_type lu` solver choice may be helpful:

    $ ./nc2petsc.py --refine 2 grn.nc grn2p5km.dat
    $ mkdir test2p5km/
    $ mpiexec -n 6 ../mahaffy -mah_read grn2p5km.dat -mah_dump test2p5km/ -snes_monitor -pc_type asm -sub_pc_type lu -snes_max_it 2000

Around 1 km grid (e.g. `--refine 1`) my 16Gb desktop runs out of memory.

To run a lower resolution, and watch the residual, do (for example)

    $ ncks -d x1,,,4 -d y1,,,4 grn.nc grn20km.nc
    $ ./nc2petsc.py grn20km.nc grn20km.dat
    $ mkdir test20km/
    $ ../mahaffy -mah_read grn20km.dat -snes_monitor -snes_vi_monitor_residual -mah_dump test20km/

