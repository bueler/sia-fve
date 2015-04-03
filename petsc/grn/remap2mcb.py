#!/usr/bin/env python
#
# (C) 2015 Ed Bueler
#
# FIXME will attempt to remap fields from pism_Greenland_5km_v1.1.nc (or
# similar) to grid used in MCdataset-2014-11-19.nc (or similar)

import argparse
import sys
import numpy as np

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
# positional
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
    nctarg = NC(args.outname, 'r')
except:
    print "ERROR: can't read from file %s ..." % args.outname
    sys.exit(12)

# read from source file
x1src = ncsrc.variables['x1'][:]
y1src = ncsrc.variables['y1'][:]
cmbsrc = np.squeeze(ncsrc.variables['climatic_mass_balance'][:])

print ncsrc.proj4
print cmbsrc.shape

# read from target file
xtarg = nctarg.variables['x'][:]
ytarg = nctarg.variables['y'][:]
bedtarg = np.squeeze(nctarg.variables['bed'][:])

print nctarg.proj4
print bedtarg.shape

