petsc/grn/mcb
=============

The SeaRISE data for bed is not very good.  The "mass-conserving beds" from
the U Irvine group is an improvement:

M. Morlighem, E. Rignot, J. Mouginot, H. Seroussi and E. Larour,
[_Deeply incised submarine glacial valleys beneath the Greenland Ice Sheet_](http://www.nature.com/ngeo/journal/vaop/ncurrent/full/ngeo2167.html),
Nat. Geosci., 7, 418-422, 2014, doi:10.1038/ngeo2167

But the data, which is available in NetCDF from NSIDC, has neither bathymetry
nor surface mass balance, and is on a different projection.  So scripts here
deal with that.

Warning: These scripts are not finished!

Do this as preprocessing:

    $ ./getmcb.sh  # downloads 1.7Gb file from NSIDC site and generates `mcb4500.nc`
    $ cp ../pism_Greenland_5km_v1.1.nc searise5km.nc
    $ ~/pism/util/nc2cdo.py searise5km.nc
    $ ./remap2mcb.py searise5km.nc mcb4500m.nc fix4500m.nc
    $ ncview -minmax all fix4500m.nc
    $ ../grn2petsc.py fix4500m.nc fix4500m.dat

Now run:

    $ (cd ../../ && make mahaffy)
    $ mkdir test/
    $ mpiexec -n 6 ../../mahaffy -mah_read fix4500m.dat -mah_D0 20.0 -mah_showdata -draw_pause 2 -snes_monitor -pc_type asm -sub_pc_type lu -snes_max_it 200 -mah_dump test/

Generate figures:

    $ ../../../figsmahaffy.py --profile --map

