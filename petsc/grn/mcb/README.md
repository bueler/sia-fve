petsc/grn/mcb
=============

The SeaRISE data for bed elevation and ice thickness is not very good.  The
"mass-conserving beds" (MCB) from the U Irvine group is an improvement:

M. Morlighem, E. Rignot, J. Mouginot, H. Seroussi and E. Larour,
[_Deeply incised submarine glacial valleys beneath the Greenland Ice Sheet_](http://www.nature.com/ngeo/journal/vaop/ncurrent/full/ngeo2167.html),
Nat. Geosci., 7, 418-422, 2014, doi:10.1038/ngeo2167

But the data, which is available in NetCDF from NSIDC, has neither bathymetry
nor surface mass balance, and is on a different projection.  So scripts here
deal with that.

Coarse resolution example
-------------------------

The first step of preprocessing is to download a 1.7Gb file `MCdataset-2014-11-19.nc`
from NSIDC site and generate a coarser version by averaging over blocks of a
certain size.  Blocks of 30 x 30 gives 30*150m = 4500m resolution:

    $ ./getmcb.sh
    $ ./average --block 30 mcb4500.nc

The next preprocessing step, still at the NetCDF level, with a merge of SeaRISE
and MCB data:

    $ cp ../pism_Greenland_5km_v1.1.nc searise5km.nc
    $ ~/pism/util/nc2cdo.py searise5km.nc
    $ ./remap2mcb.py searise5km.nc mcb4500m.nc fix4500m.nc

View and do final conversion to PETSc binary:

    $ ncview -minmax all fix4500m.nc
    $ ../grn2petsc.py fix4500m.nc fix4500m.dat

Now run `mahaffy` on it:

    $ (cd ../../ && make mahaffy)
    $ mkdir test/
    $ mpiexec -n 6 ../../mahaffy -mah_read fix4500m.dat -mah_D0 20.0 -mah_showdata -draw_pause 2 -snes_monitor -pc_type asm -sub_pc_type lu -snes_max_it 200 -mah_dump test/

Generate figures:

    $ ../../../figsmahaffy.py --profile --map

Higher resolution
-----------------

Redo the above, starting with

    $ ./getmcb.sh 10 mcb1500m.nc

for example.  (Note that if the big file `MCdataset-2014-11-19.nc` is already
present, it will not be downloaded again.)  Continue with all remaining steps.

