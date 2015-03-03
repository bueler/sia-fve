#!/usr/bin/env python
#
# (C) 2015 Ed Bueler
#
# Generate PETSc binary format file from Greenland NetCDF file.
# See README.md to set up.

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt

try:
    from netCDF4 import Dataset as NC
except:
    print "netCDF4 (netcdf4-python) is not installed"
    sys.exit(1)

try:
    import PetscBinaryIO as pbio
except:
    print "'import PetscBinaryIO' failed"
    print "need link to petsc/bin/pythonscripts/PetscBinaryIO.py?"
    sys.exit(2)

try:
    import petsc_conf
except:
    print "'import petsc_conf.py' failed"
    print "need link to petsc/bin/pythonscripts/petsc_conf.py?"
    sys.exit(2)

parser = argparse.ArgumentParser(description='Generate PETSc binary format file from Greenland NetCDF file.')
# positional
parser.add_argument('inname', metavar='INNAME',
                    help='name of NetCDF input file with topg variable (e.g. pism_Greenland_5km_v1.1.nc)',
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
thk = np.squeeze(nc.variables['thk'][:])
topg = np.squeeze(nc.variables['topg'][:])
cmb = np.squeeze(nc.variables['cmb'][:])

# linear interpolation plus padding in 1D: c = coarse and f = fine
def lininterp(c,r):
    f = np.zeros(r * len(c))
    for j in range(len(c)-1):
        for s in range(r):
            lam = float(s)/float(r)
            f[r*j + s] = (1.0 - lam) * c[j] + lam * c[j+1]
    j = len(c)-1
    for s in range(r):
        f[r*j + s] = c[j]
    return f

# bilinear interpolation plus padding in 2D: c = coarse and f = fine
def quadinterp(c,r):
    f = np.zeros((r*np.shape(c)[0],r*np.shape(c)[1]))
    for j in range(np.shape(c)[0]-1):
        for k in range(np.shape(c)[1]-1):
            for s in range(r):
                ls = float(s)/float(r)
                for t in range(r):
                    lt = float(t)/float(r)
                    f[r*j+s,r*k+t] =   (1.0-ls) * (1.0-lt) * c[j,k] \
                                     + ls       * (1.0-lt) * c[j+1,k] \
                                     + (1.0-ls) * lt       * c[j,k+1] \
                                     + ls       * lt       * c[j+1,k+1]
    j = np.shape(c)[0]-1
    for k in range(np.shape(c)[1]-1):
        for s in range(r):
            for t in range(r):
                lt = float(t)/float(r)
                f[r*j+s,r*k+t] = (1.0-lt) * c[j,k] + lt * c[j,k+1]
    k = np.shape(c)[1]-1
    for j in range(np.shape(c)[1]-1):
        for s in range(r):
            ls = float(s)/float(r)
            for t in range(r):
                f[r*j+s,r*k+t] = (1.0-ls) * c[j,k] + ls * c[j+1,k]
    j = np.shape(c)[0]-1
    k = np.shape(c)[1]-1
    for s in range(r):
        for t in range(r):
            f[r*j+s,r*k+t] = c[j,k]
    return f

if (refine != 1):
    print "refining grid to %.2f km ..." % (5.0/float(refine))
    x = lininterp(x,refine)
    y = lininterp(y,refine)
    thk = quadinterp(thk,refine)
    topg = quadinterp(topg,refine)
    cmb = quadinterp(cmb,refine)

# flatten
thk  = thk.flatten()
topg = topg.flatten()
cmb  = cmb.flatten()

print "variable lengths:  x = %d,  y = %d" % (np.shape(x)[0],np.shape(y)[0])
print "                   topg (flattened) = %d" % (np.shape(topg)[0])
print "                   cmb (flattened) = %d" % (np.shape(cmb)[0])

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

# convert to PETSc-type vecs
xvec = x.view(pbio.Vec)
yvec = y.view(pbio.Vec)
topgvec = topg.view(pbio.Vec)
cmbvec = cmb.view(pbio.Vec)

# open petsc binary file
io = pbio.PetscBinaryIO()

# write fields **in particular order**; names do not matter
io.writeBinaryFile(args.outname, [xvec,yvec,topgvec,cmbvec,])

