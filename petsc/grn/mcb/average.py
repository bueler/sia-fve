#!/usr/bin/env python
#
# (C) 2015 Ed Bueler

import argparse
import sys
import numpy as np

try:
    from netCDF4 import Dataset as NC
except:
    print "netCDF4 (=netcdf4-python) is not installed? ... see http://unidata.github.io/netcdf4-python/"
    sys.exit(1)

parser = argparse.ArgumentParser(description='Read and average 150 m NetCDF file MCdataset-2014-11-19.nc onto subgrid.')
parser.add_argument('outname', metavar='OUTNAME',
                    help='output NetCDF file to create (e.g. fix4500m.nc)',
                    default='')
parser.add_argument('--block', action='store', metavar='N', type=int,
                    help='average over N x N blocks',
                    default=30)
args = parser.parse_args()

N = args.block

inname = 'MCdataset-2014-11-19.nc'

print "reading from %s ..." % inname
try:
    nc = NC(inname, 'r')
except:
    print "ERROR: can't open %s for reading ..." % inname
    sys.exit(11)

try:
    ncout = NC(args.outname, 'w', format='NETCDF3_CLASSIC')
except:
    print "ERROR: can't open %s for writing ..." % args.outname
    sys.exit(13)

x = nc.variables['x'][:].astype(np.float64)
y = nc.variables['y'][:].astype(np.float64)
dx = x[1] - x[0]
dy = y[1] - y[0]
print "    ... which has dimensions (y,x)=(%d,%d) and dy=%.2f,dx=%.2f ..." % \
      (len(y),len(x),dy,dx)

Mx = len(x) / N
My = len(y) / N
dxout = dx * N
dyout = dy * N
print "new grid in %s has dimensions (y,x)=(%d,%d) and dy=%.2f,dx=%.2f ..." % \
      (args.outname,My,Mx,dyout,dxout)
avbed = np.zeros((My,Mx))
avthk = np.zeros((My,Mx))
xoutstart = np.average(x[0:N])
youtstart = np.average(y[0:N])
xout = np.linspace(xoutstart,xoutstart + (Mx-1) * dxout,Mx)
yout = np.linspace(youtstart,youtstart + (My-1) * dyout,My)

print "averaging over %d x %d blocks:" % (N,N)
for var in [(avbed,'bed'), (avthk,'thickness')]:
    vout = var[0]
    varname = var[1]
    if varname == 'bed':
        v = np.squeeze(nc.variables[varname][:].filled())
    else:
        v = np.squeeze(nc.variables[varname][:])
    print "    '%s', which has dimensions (y,x)=(%d,%d)" % (varname,v.shape[0],v.shape[1])
    for k in range(My):
        for j in range(Mx):
            block = v[N*k:N*(k+1),N*j:N*(j+1)].astype(np.float64)
            if np.all(block < -9000.0):
                vout[k,j] = -9999.0
            else:
                vout[k,j] = np.average(block[block>=-9000.0])

nc.close()

def defvar(nc, name, units, stdname):
    var = nc.createVariable(name, 'f4', dimensions=("y1", "x1"))
    var.units = units
    var.standard_name = stdname
    return var

print "creating x1,y1 dimensions and variables in %s ..." % args.outname
ncout.createDimension("x1", size=len(xout))
ncout.createDimension("y1", size=len(yout))
x_var = ncout.createVariable("x1", 'f4', dimensions=("x1",))
x_var.units = "m"
x_var.long_name = "easting"
x_var.standard_name = "projection_x_coordinate"
x_var[:] = xout
y_var = ncout.createVariable("y1", 'f4', dimensions=("y1",))
y_var.units = "m"
y_var.long_name = "northing"
y_var.standard_name = "projection_y_coordinate"
y_var[:] = yout

print "writing averaged thickness variable into %s ..." % args.outname
thk_var = defvar(ncout, "thk", "m", "land_ice_thickness")
thk_var[:] = avthk

print "writing averaged bed variable into %s ..." % args.outname
bed_var = defvar(ncout, "topg_nobathy", "m", "bedrock_altitude")
bed_var[:] = avbed

ncout.close()

