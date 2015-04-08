petsc/grn/
==========

The main script `grn2petsc.py` here reads SeaRISE data from NetCDF format and
converts it to PETSc binary format.  It requires the
[netcdf4-python module](https://github.com/Unidata/netcdf4-python).
The code `mahaffy` can read the specially-formatted binary file it produces.

Stage 1
-------

Get the data as slightly-cleaned NetCDF.  This stage requires
[NCO](http://nco.sourceforge.net/) and a [PISM](http://www.pism-docs.org)
preprocess script.

The dimensions in the `.nc` file must be multiples of 3
because `mahaffy.c` uses a FD-by-coloring-evaluated Jacobian on a periodic
grid with stencil width 1.  Thus we trim the x dimension from 301 values to 300,
but this removes only ocean cells.  The y dimension is already 561, which is
divisible by 3.  Also note 1 / (910 * 31556926) = 3.48228182586954e-11
is the conversion factor to turn  kg m-2 a-1  into  m s-1  for ice with
density 910 kg m-3.

Do:

    $ cd ~/pism/examples/std-greenland/
    $ ./preprocess.sh                       # downloads Greenland_5km_v1.1.nc with wget if not present
    $ cd ~/sia-fve/petsc/grn/               # back to this dir
    $ ln -s ~/pism/examples/std-greenland/pism_Greenland_5km_v1.1.nc
    $ ./cleangrn.sh

The last command uses NCO to clean up pism_Greenland_5km_v1.1.nc and creates `grn.nc`.

Stage 2
-------

Convert `grn.nc` to PETSc binary `grn5km.dat` using PETSc scripts:

    $ ln -s ~/petsc/bin/petsc-pythonscripts/PetscBinaryIO.py    # or similar for maint/3.5
    $ ln -s ~/petsc/bin/petsc-pythonscripts/petsc_conf.py
    $ python grn2petsc.py grn.nc grn5km.dat    # defaults to 5km

Note `grn2petsc.py` uses a module `interpad.py`.

Stage 3
-------

Build the executable `mahaffy`, and run on 5km grid using an over-estimate of
diffusivity:

    $ (cd ../ && make mahaffy)
    $ mkdir test/
    $ mpiexec -n 6 ../mahaffy -mah_read grn5km.dat -mah_showdata -draw_pause 2 -snes_monitor -mah_divergetryagain -mah_dump test/

This run only takes a few minutes and uses the data as is.

Generating figures:
-------------------

Generate `.pdf` and `.png` figures:

    $ cd test/
    $ ../../figsmahaffy.py --observed

Higher/lower resolution
-----------------------

A higher-res (2.5 km), parallel version might look like this; the robust
`-pc_type asm -sub_pc_type lu` solver choice may be helpful:

    $ python grn2petsc.py --refine 2 grn.nc grn2p5km.dat
    $ mkdir test2p5km/
    $ mpiexec -n 6 ../mahaffy -mah_read grn2p5km.dat -mah_dump test2p5km/ -snes_monitor -pc_type asm -sub_pc_type lu -snes_max_it 2000

Around 700m grid (e.g. `--refine 7`) my 16Gb desktop runs out of memory.

To run a lower resolution, and watch the residual, do (for example)

    $ ncks -d x1,,,4 -d y1,,,4 grn.nc grn20km.nc
    $ ./grn2petsc.py grn20km.nc grn20km.dat
    $ mkdir test20km/
    $ ../mahaffy -mah_read grn20km.dat -snes_monitor -snes_vi_monitor_residual -mah_dump test20km/

