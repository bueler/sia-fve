petsc/grn/
==========

Theses instructions document how to get
[SeaRISE](http://websrv.cs.umt.edu/isis/index.php/SeaRISE_Assessment)
data in
[NetCDF](http://www.unidata.ucar.edu/software/netcdf/)
format, do different amounts of pre-processing on it, convert it to
[PETSc binary format](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Viewer/PetscViewerBinaryOpen.html),
and run the code `mahaffy` on it.

These steps require the
[netcdf4-python](https://github.com/Unidata/netcdf4-python)
module and
[NCO](http://nco.sourceforge.net/),
in addition to stuff needed in `../`.


Get the data, link scripts, and compile
---------------------------------------

This stage requires a copy of [PISM](http://www.pism-docs.org).  In the PISM
`examples/std-greenland/` directory is a script to download the SeaRISE dataset
`Greenland_5km_v1.1.nc` and preprocess it.  Do:

    $ cd ~/pism/examples/std-greenland/
    $ ./preprocess.sh                   # downloads with wget if needed
    $ cd ~/sia-fve/petsc/grn/           # back to this dir
    $ ln -s ~/pism/examples/std-greenland/pism_Greenland_5km_v1.1.nc

Of course, we need an executable:

    $ (cd ../ && make mahaffy)


Quick start
-----------

All non-optional actions below are executed by a single bash script:

    $ ./quickstart.sh


Clean the data
--------------

The next stages all work on the data as NetCDF files, so their results can be
viewed with `ncview` or similar.

For convenience we make a local copy `grn.nc` that is smaller because it
has only the variables we want.

    $ ncks -O -v x1,y1,thk,topg,climatic_mass_balance pism_Greenland_5km_v1.1.nc grn.nc

We convert the climatic mass balance variable from units  kg m-2 a-1  to  m s-1:

    $ ./inplace.py --fixcmbunits --ranges grn.nc

We change the ocean topography to not be so deep, because the ice flow
model does not know about flotation.  And we change the ocean value of climatic
mass balance to be very negative; this simulates strong calving:

    $ ./inplace.py --oceanfix --ranges grn.nc


Optional preprocessing
----------------------

  * _Refinement to higher resolution._  To refine the grid from the native 5 km
  to higher resolution, by bilinear interpolation, we use `refine.py` which
  calls methods from `q1ops.py`.  Do:

      $ ./refine.py --factor 2 grn.nc grn2p5km.nc
      $ ./inplace.py --ranges grn2p5km.nc          # same as for grn.nc

  * _Decimation to lower resolution._  To run a lower resolution, you can use
  NCO.  For example:

      $ ncks -d x1,,,4 -d y1,,,4 grn.nc grn20km.nc
      $ ./inplace.py --ranges grn20km.nc          # different from grn.nc

  * _Apply sweeps to remove bed dips._  This seems to help just slightly, but is
  immoral:

      $ ./inplace.py --bedundip --sweeps 2 --ranges grn.nc

  * _Grid trim as needed to use FD-by-coloring Jacobian._  To use a
  finite-difference Jacobian, which PETSc calculates using graph "coloring",
  we must make the dimensions divisible by 3.  This is because we have a
  periodic grid and stencil width 1 (i.e. the stencil is 3 points across).
  Since the data includes a number of ocean cells at the
  edge of the computational domain, we can do this just by trimming the `x1`
  and/or `y1` dimensions to the next smallest multiple of 3.  First do this
  to examine the dimensions:

      $ ncdump -h grn.nc | head

  In this example we see that the `y1` dimension of 561 is already divisible
  by 3.  But we trim the `x1` dimension:

      $ ncks -O -d x1,0,299 grn.nc grn.nc


Convert to PETSC binary
-----------------------

This stage uses the links to PETSc scripts setup in `../README.md`:

    $ ../nc2petsc.py grn.nc grn.dat


Run the model
-------------

Run on the 5km grid; the `-mah_showdata` option means that the data is viewed
with PETSc's X viewer:

    $ mkdir test/
    $ mpiexec -n 6 ../mahaffy -mah_read grn.dat -cs_D0 10.0 -snes_monitor -mah_showdata -draw_pause 2 -pc_type asm -sub_pc_type lu -mah_dump test/

This run only takes a few minutes.


More optional steps
-------------------

  * _Do "manual" time steps after divergence, to approach steady state._  If we
  do

        $ ./quickstart.sh "-mah_notry"

  then we get a run that diverges at step 11 of the continuation scheme, i.e.
  it stops with a state that converges with eps=4.64e-4.  Thus it is not at
  steady state with the full SIA model.  Now do 1000 a of time steps to head
  toward steady state:

        $ mkdir testcont/
        $ mpiexec -n 6 ../mahaffy -mah_read test/unnamed.dat -mah_readinitial -cs_start 9 -mah_steps 100 -mah_dtBE 10.0 -mah_notry -mah_dump testcont/

  * _Generate figures._  Generate `.pdf` and `.png` figures:

        $ cd test/
        $ ../../figsmahaffy.py --observed

  * _Higher resolution._  Try

        $ QSSHOW= QSREFINE=2 ./quickstart.sh "-snes_max_it 200"

  Around 1 km grid (e.g. `QSREFINE=5`) my 16Gb desktop runs out of memory.

