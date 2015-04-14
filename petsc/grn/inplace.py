#!/usr/bin/env python
# (C) 2015 Ed Bueler
#
# See README.md for example usage.

import argparse
import sys
import numpy as np

from q1ops import locquadinterp

try:
    from netCDF4 import Dataset as NC
except:
    print "netCDF4 (netcdf4-python) is not installed"
    sys.exit(1)

parser = argparse.ArgumentParser(description='In-place operations on .nc with ice sheet data.')
parser.add_argument('name', metavar='NAME',
                    help='NetCDF file with x1,y1,topg,cmb,thk variables (e.g. grn.nc); on output has has topg-->topg_orig and new topg')
parser.add_argument('--bedundip', action='store_true',
                    help='sweeps which damp-out concave-up places on bed')
parser.add_argument('--fixcmbunits', action='store_true',
                    help='convert cmb from  kg m-2 year-1  to  m s-1')
parser.add_argument('--oceanfix', action='store_true',
                    help='modify cmb and topg in ocean')
parser.add_argument('--ranges', action='store_true',
                    help='print out dimensions and variable ranges')
parser.add_argument('--sweeps', action='store', metavar='N', default=1, type=int,
                    help='number of "--bedundip" sweeps; default=%(default)d)')
args = parser.parse_args()

try:
    nc = NC(args.name, 'a')
except:
    print "ERROR: can't open file %s for reading and writing ..." % args.name
    sys.exit(11)

def defvar(nc, name, units, stdname):
    var = nc.createVariable(name, 'f4', dimensions=("y1", "x1"))
    var.units = units
    var.standard_name = stdname
    return var

secpera = 31556926.0

if args.fixcmbunits:
    print "creating new cmb variable with units  m s-1  ..."
    cmb_var = defvar(nc, "cmb", "m s-1", "")
    factor = 1.0 / (910.0 * secpera) #= 3.48228182586954e-11
    cmb_var[:] = factor * nc.variables['climatic_mass_balance'][:]

if args.oceanfix:
    print "applying ocean fixes ..."
    topg = np.squeeze(nc.variables['topg'])
    cmb  = np.squeeze(nc.variables['cmb'])
    thk  = np.squeeze(nc.variables['thk'])
    for b,m,H in np.nditer([topg,cmb,thk], op_flags=['readwrite']):
         if (H <= 0.0):
             if (b < -250.0):
                 b[...] = -250.0
             if (b < -200.0):
                 m[...] = -30.0 / 31556926.0
             elif (b < -100.0):
                 m[...] = -10.0 / 31556926.0
             elif (b < -50.0):
                 m[...] = -5.0 / 31556926.0
    print "overwriting topg and cmb variable with modified versions ..."
    nc.variables['topg'][:] = topg[:]
    nc.variables['cmb'][:] = cmb[:]

damp = 0.5
def dampme(f, v):
    av = locquadinterp(f,0.5,0.5)
    if v < av:
        return v + damp * (av - v)
    else:
        return v

if args.bedundip:
    print "applying sweeps to dampen bed concave-up areas ..."
    topg = np.squeeze(nc.variables['topg'])
    for N in range(args.sweeps):
        newtopg = topg.copy()
        for j in range(0,np.shape(topg)[0]-3,2):
            for k in range(0,np.shape(topg)[1]-3,2):
                newtopg[j+1][k+1] = dampme(topg[j:j+3:2,  k:k+3:2],  topg[j+1][k+1])
                newtopg[j+1][k+2] = dampme(topg[j:j+3:2,  k+1:k+4:2],topg[j+1][k+2])
                newtopg[j+2][k+1] = dampme(topg[j+1:j+4:2,k:k+3:2],  topg[j+2][k+1])
                newtopg[j+2][k+2] = dampme(topg[j+1:j+4:2,k+1:k+4:2],topg[j+2][k+2])
        topg = newtopg.copy()
        print ".",
        sys.stdout.flush()
    print ""
    print "overwriting topg with modified version ..."
    nc.variables['topg'][:] = topg[:]

if args.ranges:
    x = nc.variables['x1'][:]
    y = nc.variables['y1'][:]
    topg = np.squeeze(nc.variables['topg'])
    cmb  = np.squeeze(nc.variables['cmb'])
    thk  = np.squeeze(nc.variables['thk'])
    print "dimensions:  x1 = %d,  y1 = %d;  variable ranges:" % (np.shape(x)[0],np.shape(y)[0])
    print "      %10.4f <= topg <= %10.4f  (m)" % (topg.min(), topg.max())
    print "      %10.3e <= cmb  <= %10.3e  (m a-1)" % (cmb.min()*secpera, cmb.max()*secpera)
    print "      %10.4f <= thk  <= %10.4f  (m)" % (thk.min(), thk.max())

nc.close()

