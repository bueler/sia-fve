#!/usr/bin/env python
# (C) 2015 Ed Bueler
#
# See README.md for example usage.

import argparse
import sys
import numpy as np

from q1ops import lininterp, quadinterp

try:
    from netCDF4 import Dataset as NC
except:
    print "netCDF4 (netcdf4-python) is not installed"
    sys.exit(1)

parser = argparse.ArgumentParser(description='Refine a grid by a given factor.')
parser.add_argument('inname', metavar='INNAME',
                    help='input NetCDF file with x1,y1,topg,cmb,thk variables (e.g. grn.nc)')
parser.add_argument('outname', metavar='OUTNAME',
                    help='output NetCDF file with refined grid and same variables')
parser.add_argument('--factor', action='store', metavar='N', default=2, type=int,
                    help='factor by which grid is refined; default=%(default)d)')
args = parser.parse_args()

try:
    nc = NC(args.inname, 'r')
except:
    print "ERROR: can't open %s for reading ..." % args.inname
    sys.exit(11)

x1in = nc.variables['x1']
y1in = nc.variables['y1']
dx = x1in[1] - x1in[0]
dy = y1in[1] - y1in[0]
print "original grid in %s has spacing %.3f km in x and %.3f km in y ..." % \
      (args.inname,dx/1000.0,dy/1000.0)

print "refining vars x1,y1,topg,cmb,thk to %.3f km ..." % ((dx/1000.0)/float(args.factor))
x1 = lininterp(x1in,args.factor)
y1 = lininterp(y1in,args.factor)
topg = quadinterp(np.squeeze(nc.variables['topg']),args.factor)
cmb = quadinterp(np.squeeze(nc.variables['cmb']),args.factor)
thk = quadinterp(np.squeeze(nc.variables['thk']),args.factor)

nc.close()

try:
    nc = NC(args.outname, 'w')
except:
    print "ERROR: can't open %s for writing ..." % args.outname
    sys.exit(12)

def defvar(nc, name, units, stdname):
    var = nc.createVariable(name, 'f4', dimensions=("y1", "x1"))
    var.units = units
    var.standard_name = stdname
    return var

print "writing dimensions and axes to %s ..." % args.outname
nc.createDimension("x1", size=len(x1))
nc.createDimension("y1", size=len(y1))
x1_var = nc.createVariable("x1", 'f4', dimensions=("x1",))
x1_var.units = "m"
x1_var.long_name = "easting"
x1_var.standard_name = "projection_x_coordinate"
x1_var[:] = x1
y1_var = nc.createVariable("y1", 'f4', dimensions=("y1",))
y1_var.units = "m"
y1_var.long_name = "northing"
y1_var.standard_name = "projection_y_coordinate"
y1_var[:] = y1

print "writing variables to %s ..." % args.outname
topg_var = defvar(nc, "topg", "m", "bedrock_altitude")
topg_var[:] = topg
cmb_var = defvar(nc, "cmb", "m s-1", "")
cmb_var[:] = cmb
thk_var = defvar(nc, "thk", "m", "land_ice_thickness")
thk_var[:] = thk

nc.close()
print "done."

