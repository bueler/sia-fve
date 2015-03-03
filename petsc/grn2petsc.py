#!/usr/bin/env python
#
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
# positional
parser.add_argument('inname', metavar='INNAME',
                    help='name of NetCDF input file with x1,y1,topg,cmb,thk variables (e.g. pism_Greenland_5km_v1.1.nc)',
                    default='')
parser.add_argument('outname', metavar='OUTNAME',
                    help='name of output PETSc binary file (e.g. grn.dat)',
                    default='')
# optional
parser.add_argument('--refine', action='store', metavar='N',
                    help='refine by factor of N (N=1 is no change, N=2 for 2.5km, N=5 for 5km)',
                    default=1)
args = parser.parse_args()
refine = int(args.refine)

try:
    nc = NC(args.inname, 'r')
except:
    print "ERROR: can't read from file %s ..." % args.inname
    sys.exit(11)

# read axes
x = nc.variables['x1'][:]
y = nc.variables['y1'][:]

# load data
topg = np.squeeze(nc.variables['topg'][:])
cmb = np.squeeze(nc.variables['cmb'][:])
thk = np.squeeze(nc.variables['thk'][:])

if (refine != 1):
    print "refining grid to %.2f km ..." % (5.0/float(refine))
    x = lininterp(x,refine)
    y = lininterp(y,refine)
    topg = quadinterp(topg,refine)
    cmb = quadinterp(cmb,refine)
    thk = quadinterp(thk,refine)

# flatten
topg = topg.flatten()
cmb  = cmb.flatten()
thk  = thk.flatten()

print "variable lengths:  x = %d,  y = %d" % (np.shape(x)[0],np.shape(y)[0])
print "                   topg (flattened) = %d" % (np.shape(topg)[0])
print "                   cmb (flattened) = %d" % (np.shape(cmb)[0])
print "                   thk (flattened) = %d" % (np.shape(thk)[0])

def showranges(t,c):
    print "      %10.4f <= topg <= %10.4f  (m)" % (t.min(), t.max())
    print "      %10.3e <= cmb  <= %10.3e  (m s-1)" % (c.min(), c.max())

# modify b and cmb in ocean
print "before ocean fixes:"
showranges(topg,cmb)
for j in range(len(cmb)):
     if (thk[j] <= 0.0):
         if (topg[j] < -250.0):
             topg[j] = -250.0
         if (topg[j] < -200.0):
             cmb[j] = -30.0 / 31556926.0
         elif (topg[j] < -100.0):
             cmb[j] = -10.0 / 31556926.0
         elif (topg[j] < -50.0):
             cmb[j] = -5.0 / 31556926.0
print "after ocean fixes:"
showranges(topg,cmb)

print "observed thickness:"
print "      %10.4f <= thk  <= %10.4f  (m)" % (thk.min(), thk.max())

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

