#!/usr/bin/env python
# (C) 2015 Ed Bueler

import argparse
import sys
import numpy as np

from interpad import quadinterp

try:
    from netCDF4 import Dataset as NC
except:
    print "netCDF4 (netcdf4-python) is not installed"
    sys.exit(1)

parser = argparse.ArgumentParser(description='A special kind of bed-smoothing sweep: damp-out concave-up places on bed.')
parser.add_argument('name', metavar='NAME',
                    help='NetCDF file with x1,y1,topg,cmb,thk variables (e.g. grn.nc); on output has has topg-->topg_orig and new topg')
parser.add_argument('--sweeps', action='store', metavar='N', default=1, type=int,
                    help='number of sweeps; default=%(default)d)')
args = parser.parse_args()

try:
    nc = NC(args.name, 'a')
except:
    print "ERROR: can't open file %s for reading and writing ..." % args.inname
    sys.exit(11)

print "reading axes x,y and fields topg,cmb,thk ..."
x = nc.variables['x1'][:]
y = nc.variables['y1'][:]
topg = np.squeeze(nc.variables['topg'][:])
cmb = np.squeeze(nc.variables['cmb'][:])
thk = np.squeeze(nc.variables['thk'][:])

print "dimensions:  x = %d,  y = %d" % (np.shape(x)[0],np.shape(y)[0])

def showranges():
    secpera = 31556926.0
    print "      %10.4f <= topg <= %10.4f  (m)" % (topg.min(), topg.max())
    print "      %10.3e <= cmb  <= %10.3e  (m a-1)" % (cmb.min()*secpera, cmb.max()*secpera)
    print "      %10.4f <= thk  <= %10.4f  (m)" % (thk.min(), thk.max())

print "variable ranges:"
showranges()

def defvar(nc, name, units, stdname):
    var = nc.createVariable(name, 'f4', dimensions=("y1", "x1"))
    var.units = units
    var.standard_name = stdname
    return var

print "creating topg_orig variable and copying  topg --> topg_orig  ..."
orig_var = defvar(nc, "topg_orig", "m", "")
orig_var[:] = topg

print "overwriting topg variable with modified ..."
nc.variables['topg'][:] = 1000.0 * topg[:]  # FIXME: test version

nc.close()

