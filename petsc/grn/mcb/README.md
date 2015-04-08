petsc/grn/mcb
=============

The SeaRISE data for bed elevation and ice thickness is not very good.  The
"mass-conserving bed" (MCB) from the U Irvine group is an improvement:

M. Morlighem, E. Rignot, J. Mouginot, H. Seroussi and E. Larour,
[_Deeply incised submarine glacial valleys beneath the Greenland Ice Sheet_](http://www.nature.com/ngeo/journal/vaop/ncurrent/full/ngeo2167.html),
Nat. Geosci., 7, 418-422, 2014, doi:10.1038/ngeo2167

This bed data, which is available in NetCDF from NSIDC, does not have
bathymetry.  There is no surface mass balance on this same projection, that I
know of.  So scripts here deal with those issues.

Coarse resolution example
-------------------------

The first step of preprocessing is to download a 1.7Gb file `MCdataset-2014-11-19.nc`
from NSIDC site and generate a coarser version by averaging over blocks of a
certain size.  Blocks of 30 x 30 gives 30*150m = 4500m resolution:

    $ ./getmcb.sh
    $ ./average.py --block 30 mcb4500m.nc

The next preprocessing step, still at the NetCDF level, is a merge of SeaRISE
and MCB data.  In this script `remap2mcb.py` we need to hand the conversion
between projections:

    $ cp ../pism_Greenland_5km_v1.1.nc searise5km.nc
    $ ~/pism/util/nc2cdo.py searise5km.nc
    $ ./remap2mcb.py searise5km.nc mcb4500m.nc fix4500m.nc

Finally we do a conversion to PETSc binary.  See `../README.md` to set up links
so that `../grn2petsc.py` will work:

    $ ../grn2petsc.py fix4500m.nc fix4500m.dat

Now run `mahaffy` on it; these defaults do well:

    $ (cd ../../ && make mahaffy)   # if needed
    $ mkdir test/
    $ mpiexec -n 6 ../../mahaffy -mah_read fix4500m.dat -mah_showdata -draw_pause 2 -snes_monitor -pc_type asm -sub_pc_type lu -snes_max_it 200 -mah_dump test/

Generate figures on the results in `test/` as described in `../README.md`.

Higher resolution
-----------------

Redo the above, for example:

    $ ./getmcb.sh
    $ ./average.py --block 10 mcb1500m.nc
    $ ./remap2mcb.py searise5km.nc mcb1500m.nc fix1500m.nc
    $ ../grn2petsc.py fix1500m.nc fix1500m.dat
    $ mkdir test1500m/
    $ mpiexec -n 6 ../../mahaffy -mah_read fix1500m.dat -snes_monitor -pc_type asm -sub_pc_type lu -snes_max_it 200 -mah_dump test1500m/

Note that if the big file `MCdataset-2014-11-19.nc` is already present,
re-running `getmcb.sh` will not cause it to be downloaded again.

