#!/usr/bin/env python
# (C) 2015 Ed Bueler

import argparse
import sys
import numpy as np

try:
    from netCDF4 import Dataset as NC
except:
    print "netCDF4 (netcdf4-python) is not installed"
    sys.exit(1)

try:
    import PetscBinaryIO as pbio
except:
    print "'import PetscBinaryIO' failed"
    print "need link to petsc/bin/petsc-pythonscripts/PetscBinaryIO.py?"
    sys.exit(2)

try:
    import petsc_conf
except:
    print "'import petsc_conf.py' failed"
    print "need link to petsc/bin/petsc-pythonscripts/petsc_conf.py?"
    sys.exit(2)

parser = argparse.ArgumentParser(description='Generate NetCDF file from PETSc binary file.  Files have ice-sheet specific variables.')
parser.add_argument('inname', metavar='INNAME',
                    help='input PETSc binary file with x1,y1,topg,cmb,thk variables (e.g. grn.dat)')
parser.add_argument('outname', metavar='OUTNAME',
                    help='output NetCDF file (e.g. grn.nc)')
args = parser.parse_args()

try:
    nc = NC(args.outname, 'w')
except:
    print "ERROR: can't open NetCDF file %s for writing ..." % args.outname
    sys.exit(11)

def typeerror(name):
    print "ERROR: unexpected objecttype when reading '%s' ..." % name
    sys.exit(71)

print "reading axes x1,y1 from %s ..." % args.inname

io = pbio.PetscBinaryIO()
fh = open(args.inname)
x1type = io.readObjectType(fh)
if x1type == 'Vec':
    x1 = io.readVec(fh)
else:
    typeerror('x1')
y1type = io.readObjectType(fh)
if y1type == 'Vec':
    y1 = io.readVec(fh)
else:
    typeerror('y1')

Nx = len(x1)
Ny = len(y1)

def readfield(name):
    vtype = io.readObjectType(fh)
    if vtype == 'Vec':
        v = np.reshape(io.readVec(fh),(Ny,Nx))
    else:
        typeerror(name)
    return v

print "reading fields topg,cmb,thk from %s ..." % args.inname
topg = readfield('topg')
cmb  = readfield('cmb')
thk  = readfield('thk')

print "writing vars x1,y1,topg,cmb,thk into %s ..." % args.outname

def defvar(nc, name, units, stdname):
    var = nc.createVariable(name, 'f4', dimensions=("y1", "x1"))
    var.units = units
    var.standard_name = stdname
    return var

print "writing dimensions and axes to %s ..." % args.outname
nc.createDimension("x1", size=Nx)
nc.createDimension("y1", size=Ny)
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
