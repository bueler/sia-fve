#!/usr/bin/env python
# (C) 2015 Ed Bueler
#
# Generate PETSc binary format file from Greenland NetCDF file.
# See README.md to set up.

import argparse
import sys
import numpy as np

from interpad import lininterp, quadinterp

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

parser = argparse.ArgumentParser(description='Generate PETSc binary format file from Greenland NetCDF file.')
parser.add_argument('inname', metavar='INNAME',
                    help='input NetCDF file with x1,y1,topg,cmb,thk variables (e.g. grn.nc)')
parser.add_argument('outname', metavar='OUTNAME',
                    help='output PETSc binary file (e.g. grn1km.dat if --refine 5)')
parser.add_argument('--refine', action='store', metavar='N', default=1, type=int,
                    help='refine by factor of N (N=1 is no change, N=2 for 2.5km, N=5 for 1km; default=%(default)d)')
parser.add_argument('--nooceanfix', action='store_true',
                    help='DO NOT modify cmb and topg in ocean; note ocean fixes are applied AFTER refinement')
fixes')
args = parser.parse_args()
refine = int(args.refine)

try:
    nc = NC(args.inname, 'r')
except:
    print "ERROR: can't read from file %s ..." % args.inname
    sys.exit(11)

print "reading axes x,y and fields topg,cmb,thk ..."
x = nc.variables['x1'][:]
y = nc.variables['y1'][:]
topg = np.squeeze(nc.variables['topg'][:])
cmb = np.squeeze(nc.variables['cmb'][:])
thk = np.squeeze(nc.variables['thk'][:])

if (refine != 1):
    print "refining grid to %.3f km ..." % (5.0/float(refine))
    x = lininterp(x,refine)
    y = lininterp(y,refine)
    topg = quadinterp(topg,refine)
    cmb = quadinterp(cmb,refine)
    thk = quadinterp(thk,refine)

print "dimensions:  x = %d,  y = %d" % (np.shape(x)[0],np.shape(y)[0])

def showranges():
    secpera = 31556926.0
    print "      %10.4f <= topg <= %10.4f  (m)" % (topg.min(), topg.max())
    print "      %10.3e <= cmb  <= %10.3e  (m a-1)" % (cmb.min()*secpera, cmb.max()*secpera)
    print "      %10.4f <= thk  <= %10.4f  (m)" % (thk.min(), thk.max())

print "variable ranges:"
showranges()

if not args.nooceanfix:
    # modify topg and cmb in ocean
    print "applying ocean fixes ..."
    for H,b,m in np.nditer([thk,topg,cmb], op_flags=['readwrite']):
         if (H <= 0.0):
             if (b < -250.0):
                 b[...] = -250.0
             if (b < -200.0):
                 m[...] = -30.0 / 31556926.0
             elif (b < -100.0):
                 m[...] = -10.0 / 31556926.0
             elif (b < -50.0):
                 m[...] = -5.0 / 31556926.0
    print "new variable ranges:"
    showranges()

# flatten
topg = topg.flatten()
cmb  = cmb.flatten()
thk  = thk.flatten()

# convert to PETSc-type vecs
xvec = x.view(pbio.Vec)
yvec = y.view(pbio.Vec)
topgvec = topg.view(pbio.Vec)
cmbvec = cmb.view(pbio.Vec)
thkvec = thk.view(pbio.Vec)

# open petsc binary file
io = pbio.PetscBinaryIO()

# write fields **in particular order**; names do not matter
print "writing vars x,y,topg,cmb,thk into %s ..." % args.outname
io.writeBinaryFile(args.outname, [xvec,yvec,topgvec,cmbvec,thkvec,])

