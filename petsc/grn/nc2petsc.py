#!/usr/bin/env python
# (C) 2015 Ed Bueler
#
# Generate PETSc binary format file from Greenland NetCDF file.
# See README.md.

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

parser = argparse.ArgumentParser(description='Generate PETSc binary format file from Greenland NetCDF file, with optional refinement.')
parser.add_argument('inname', metavar='INNAME',
                    help='input NetCDF file with x1,y1,topg,cmb,thk variables (e.g. grn.nc)')
parser.add_argument('outname', metavar='OUTNAME',
                    help='output PETSc binary file (e.g. grn1km.dat if --refine 5)')
parser.add_argument('--refine', action='store', metavar='N', default=1, type=int,
                    help='refine by factor of N (N=1 is no change, N=2 for 2.5km, N=5 for 1km; default=%(default)d)')
args = parser.parse_args()
refine = int(args.refine)

try:
    nc = NC(args.inname, 'r')
except:
    print "ERROR: can't read from file %s ..." % args.inname
    sys.exit(11)

print "reading axes x1,y1 and fields topg,cmb,thk ..."
x1 = nc.variables['x1'][:]
y1 = nc.variables['y1'][:]
topg = np.squeeze(nc.variables['topg'][:])
cmb = np.squeeze(nc.variables['cmb'][:])
thk = np.squeeze(nc.variables['thk'][:])

nc.close()

# optionally refine
if (refine != 1):
    print "refining grid to %.3f km ..." % (5.0/float(refine))
    x1 = lininterp(x1,refine)
    y1 = lininterp(y1,refine)
    topg = quadinterp(topg,refine)
    cmb = quadinterp(cmb,refine)
    thk = quadinterp(thk,refine)

# convert to PETSc-type vecs
x1vec = x1.view(pbio.Vec)
y1vec = y1.view(pbio.Vec)
topgvec = topg.flatten().view(pbio.Vec)
cmbvec = cmb.flatten().view(pbio.Vec)
thkvec = thk.flatten().view(pbio.Vec)

# open petsc binary file
io = pbio.PetscBinaryIO()

# write fields **in a particular order**  (though names do not matter)
print "writing vars x1,y1,topg,cmb,thk into %s ..." % args.outname
io.writeBinaryFile(args.outname, [x1vec,y1vec,topgvec,cmbvec,thkvec,])

