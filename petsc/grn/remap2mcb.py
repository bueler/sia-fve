#!/usr/bin/env python
#
# (C) 2015 Ed Bueler
#
# FIXME will attempt to remap fields from pism_Greenland_5km_v1.1.nc (or
# similar) to grid used in MCdataset-2014-11-19.nc (or similar)

import argparse
import sys
import numpy as np
import matplotlib.mlab as ml
import matplotlib.pyplot as plt

try:
    from netCDF4 import Dataset as NC
except:
    print "netCDF4 (=netcdf4-python) is not installed? ... see http://unidata.github.io/netcdf4-python/"
    sys.exit(1)

try:
    from pyproj import Proj, transform
except:
    print "pyproj is not installed? ... see http://jswhit.github.io/pyproj/"
    sys.exit(2)

parser = argparse.ArgumentParser(description='FIXME')
parser.add_argument('inname', metavar='INNAME',
                    help='input file (e.g. pism_Greenland_5km_v1.1.nc)',
                    default='')
parser.add_argument('outname', metavar='OUTNAME',
                    help='output file to modify by adding remapped fields from input (e.g. mcb4500m.nc)',
                    default='')
args = parser.parse_args()

try:
    ncsrc = NC(args.inname, 'r')
except:
    print "ERROR: can't read from file %s ..." % args.inname
    sys.exit(11)

try:
    nctarg = NC(args.outname, 'a')
except:
    print "ERROR: can't read from file %s ..." % args.outname
    sys.exit(12)

# read from source file
x1src = ncsrc.variables['x1'][:]  # is increasing
y1src = ncsrc.variables['y1'][:]
y1src = y1src[::-1]  #  is now decreasing  FIXME presumably this means y dim in cmb should change too
cmbsrc = np.squeeze(ncsrc.variables['climatic_mass_balance'][:])
print "source grid from %s has dimensions (%d,%d)" % (args.inname,cmbsrc.shape[0],cmbsrc.shape[1])

# NOTE  ncsrc.proj4 is type 'unicode'
projsrc = Proj(ncsrc.proj4.encode('ascii','ignore')) # type 'str'

# read from target file
xtarg = nctarg.variables['x'][:]  # is increasing
ytarg = nctarg.variables['y'][:]  # is decreasing
bedtarg = np.squeeze(nctarg.variables['bed'][:])
print "target grid from %s has dimensions (%d,%d)" % (args.outname,bedtarg.shape[0],bedtarg.shape[1])

# NOTE  nctarg.proj4 = +init=espg:3413 is misspelled
projtarg = Proj('+init=epsg:3413')

def insrc(x,y):
    return (x >= x1src.min()) & (x <= x1src.max()) & (y >= y1src.min()) & (y <= y1src.max())

def instr(flg):
    if flg:
        return "in source grid"
    else:
        return "outside"

def demoremap():
    # show where four corners and middle and more points of target will go in src region
    ix = [  0,   0, 333, 333, 166,  50,  50, 300, 300]
    iy = [  0, 598,   0, 598, 298,  50, 550,  50, 550]
    for k in range(len(ix)):
        xout, yout = transform(projtarg, projsrc, xtarg[ix[k]], ytarg[iy[k]])
        flg = insrc(xout, yout)
        print "(%11.2f,%11.2f) --> (%11.2f,%11.2f)  ... %s" % \
              (xtarg[ix[k]], ytarg[iy[k]], xout, yout, instr(flg))
        if flg:
            iix = ml.find(x1src <= xout).max()
            iiy = ml.find(y1src <= yout).min()
            print "    box is [%11.2f,%11.2f] x [%11.2f,%11.2f]" % \
                  (x1src[iix],x1src[iix+1],y1src[iiy],y1src[iiy-1])

#demoremap()

cmbtarg = np.zeros(bedtarg.shape)
for j in range(len(xtarg)):
    for k in range(len(ytarg)):
        FIXME


#def def_var(nc, name, units, fillvalue):
#    # dimension transpose is standard: "float thk(y, x)" in NetCDF file
#    var = nc.createVariable(name, 'f', dimensions=("y", "x"), fill_value=fillvalue)
#    var.units = units
#    return var
#smb_var = def_var(nc, "climatic_mass_balance", "kg m-2 s-1", fill_value)
#smb_var.standard_name = "land_ice_surface_specific_mass_balance_flux"
#smb_var[:] = smb



ncsrc.close()
nctarg.close()
