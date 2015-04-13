petsc/grn/
==========

Theses instructions document how to read
[SeaRISE](http://websrv.cs.umt.edu/isis/index.php/SeaRISE_Assessment)
data from
[NetCDF](http://www.unidata.ucar.edu/software/netcdf/)
format, do different amounts of pre-processing on it, and convert it to
[PETSc binary format](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Viewer/PetscViewerBinaryOpen.html).

These steps require the
[netcdf4-python](https://github.com/Unidata/netcdf4-python)
module and
[NCO](http://nco.sourceforge.net/).

The code `mahaffy` can read the specially-formatted `.dat` binary file we produce.


Get the data from the PISM source
---------------------------------

This stage requires a copy of [PISM](http://www.pism-docs.org).  In the PISM
examples are scripts to download the SeaRISE dataset `Greenland_5km_v1.1.nc`
and preprocess it.  First do:

    $ cd ~/pism/examples/std-greenland/
    $ ./preprocess.sh                   # downloads with wget if needed
    $ cd ~/sia-fve/petsc/grn/           # back to this dir

For convenience we also make a local copy `grn.nc` that is smaller because it
has only the variables we want.

    $ ncks -v x1,y1,thk,topg,climatic_mass_balance ~/pism/examples/std-greenland/pism_Greenland_5km_v1.1.nc grn.nc


Clean the data
--------------

These stages all work on the data as NetCDF files, so their results can be
viewed with `ncview` or similar.  First we we convert the climatic mass balance
variable from units  kg m-2 a-1  to  m s-1:

    $ ./inplace.py --fixcmbunits --ranges grn.nc

Also we change the ocean topography to not be so deep, because the ice flow
model does not know about flotation.  And we change the ocean value of climatic
mass balance to very negative, because this simulates strong calving.

    $ ./inplace.py --oceanfix --ranges grn.nc

FIXME: not implemented:
Optionally we apply sweeps to remove bed dips:

    $ ./inplace.py --bedundip --sweeps 2 --ranges grn.nc


Optional refinement to higher resolution
----------------------------------------

To refine the grid from the native 5 km to higher resolution, by bilinear
interpolation, do:

    $ ./refine.py --factor 2 grn.nc grn2p5km.nc
    $ ./inplace.py --ranges grn2p5km.nc          # should be same as for grn.nc


Optional lower resolution
-------------------------

To run a lower resolution, you can use NCO.  For example:

    $ ncks -d x1,,,4 -d y1,,,4 grn.nc grn20km.nc
    $ ./nc2petsc.py grn20km.nc grn20km.dat
    $ mkdir test20km/
    $ ../mahaffy -mah_read grn20km.dat -snes_monitor -mah_dump test20km/


If you want to use FD-by-coloring Jacobian
------------------------------------------

_Most users can ignore these steps_

To use a finite-difference Jacobian, which PETSc calculates using graph
"coloring", care must be taken to ensure the dimensions are divisible by 3.
This is because we have a periodic grid and stencil width 1 (i.e. the stencil is
3 points across).

Since the data includes a number of ocean cells at the edge of the computational
domain, we can do this just by trimming the `x1` and/or `y1` dimensions to the
next smallest multiple of 3.  So first do this to examine the dimensions:

    $ ncdump -h grn.nc | head

In this example we see that the `y1` dimension of 561 is already divisible by 3.
But we trim the `x1` dimension:

    $ ncks -O -d x1,0,299 grn.nc grn.nc


Convert to PETSC binary
-----------------------

This stage uses PETSc scripts, with optional refinement.

    $ ln -s ~/petsc/bin/petsc-pythonscripts/PetscBinaryIO.py  # or similar for maint/3.5
    $ ln -s ~/petsc/bin/petsc-pythonscripts/petsc_conf.py
    $ ./nc2petsc.py grn.nc grn.dat


Run the model
-------------

Build the executable `mahaffy`, and run on 5km grid using an over-estimate of
diffusivity:

    $ (cd ../ && make mahaffy)
    $ mkdir test/
    $ mpiexec -n 6 ../mahaffy -mah_read grn.dat -mah_showdata -draw_pause 2 -snes_monitor -mah_dump test/

This run only takes a few minutes and uses the data as is.


Generate figures
----------------

Generate `.pdf` and `.png` figures:

    $ cd test/
    $ ../../figsmahaffy.py --observed


Solver options for higher resolution
------------------------------------

The robust `-pc_type asm -sub_pc_type lu` solver choice may be helpful:

    $ mkdir test2p5km/
    $ ./nc2petsc.py grn2p5km.nc grn2p5km.dat
    $ mpiexec -n 6 ../mahaffy -mah_read grn2p5km.dat -mah_dump test2p5km/ -snes_monitor -pc_type asm -sub_pc_type lu -snes_max_it 2000

Around 1 km grid (e.g. `--refine 1`) my 16Gb desktop runs out of memory.

