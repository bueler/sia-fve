petsc/grn/mcb
=============

The SeaRISE data used in `../` for bed elevation and ice thickness is not very
good.  The "mass-conserving bed" (MCB) from the U Irvine group is an improvement:

M. Morlighem, E. Rignot, J. Mouginot, H. Seroussi and E. Larour,
[_Deeply incised submarine glacial valleys beneath the Greenland Ice Sheet_](http://www.nature.com/ngeo/journal/vaop/ncurrent/full/ngeo2167.html),
Nat. Geosci., 7, 418-422, 2014, doi:10.1038/ngeo2167

This bed data, which is available in NetCDF from
[NSIDC](https://nsidc.org/data)
at the specific address given in `getmcb.sh`, does not have bathymetry.  Also
there is no surface mass balance on this same projection.  Scripts here deal
with those issues.


Get large data file
-------------------

The first step is to download a 1.7Gb file `MCdataset-2014-11-19.nc` from the
NSIDC site:

    $ ./getmcb.sh    # downloads with wget

Note that if the big file `MCdataset-2014-11-19.nc` is already present,
re-running `getmcb.sh` will _not_ cause it to be downloaded again.


Preprocess SeaRISE data
-----------------------

We get the SeaRISE data on a 5km grid and preprocess it mostly as described in
`../README.md`, but adding the `mapping` variable and fixing the projection
information:

    $ ncks -O -v x1,y1,thk,topg,climatic_mass_balance,mapping ../pism_Greenland_5km_v1.1.nc searise5km.nc
    $ ../inplace.py --fixcmbunits --ranges searise5km.nc
    $ ../inplace.py --oceanfix --ranges searise5km.nc
    $ ~/pism/util/nc2cdo.py searise5km.nc


Quick start
-----------

Do the above steps.  Then this will complete the below non-optional steps:

    $ ./run.sh 30 4500m 6 "-mah_notry"


Coarsen and merge
-----------------

At this point we have not tried using the 150 m grid directly as the
computational grid!  So we generate a coarser version by averaging over blocks
of a certain size.  Blocks of 30 x 30 gives 30*150m = 4500m resolution:

    $ ./average.py --block 30 mcb4500m.nc

The next preprocessing steps, still at the NetCDF level, do a merge of SeaRISE
and MCB data.  The script `remap2mcb.py` handles the conversion between
projections:

    $ ./remap2mcb.py searise5km.nc mcb4500m.nc fix4500m.nc

Finally we do a conversion to PETSc binary:

    $ ../../nc2petsc.py fix4500m.nc fix4500m.dat


Coarse resolution run
---------------------

Now run `mahaffy` on it; these defaults do o.k.:

    $ mkdir test4500m/
    $ mpiexec -n 6 ../../mahaffy -mah_read fix4500m.dat -cs_D0 10.0 -mah_showdata -draw_pause 2 -snes_monitor -pc_type asm -sub_pc_type lu -mah_dump test4500m/ -mah_notry


Optional actions
-----------------

  * _Generate figures._

        $ cd test4500m/
        $ ../../../figsmahaffy.py

  * _Higher resolution._  For example:

        $ ./run.sh 20 3000m 6 "-mah_notry"

  or

        $ ./run.sh 10 1500m 6 "-snes_max_it 200 -mah_dtrecovery 10.0"

  * _Apply sweeps to remove bed dips._  This is immoral:

        $ cp fix4500m.nc fixbud4500m.nc
        $ ../inplace.py --bedundip --sweeps 2 --ranges fixbud4500m.nc
        $ ../nc2petsc.py fixbud4500m.nc fixbud4500m.dat

