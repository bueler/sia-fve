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

parser = argparse.ArgumentParser(description='Generate PETSc binary file from NetCDF file.  Files have ice-sheet specific variables.')
parser.add_argument('inname', metavar='INNAME',
                    help='input NetCDF file with x1,y1,topg,cmb,thk variables (e.g. grn.nc)')
parser.add_argument('outname', metavar='OUTNAME',
                    help='output PETSc binary file (e.g. grn.dat)')
parser.add_argument("--initialisobserved", action="store_true",
                    help="instead of putting zero in for initial thickness, copy observed thickness")
args = parser.parse_args()

try:
    nc = NC(args.inname, 'r')
except:
    print "ERROR: can't open NetCDF file %s for reading  ..." % args.inname
    sys.exit(11)

print "reading axes x1,y1 and fields topg,cmb,thk_obs ..."
x1 = nc.variables['x1'][:]
y1 = nc.variables['y1'][:]
topg = np.squeeze(nc.variables['topg'][:])
cmb = np.squeeze(nc.variables['cmb'][:])
thk_obs = np.squeeze(nc.variables['thk'][:])

if args.initialisobserved:
    print "setting thk_initial = thk_obs ..."
    thk_init = thk_obs.copy()
else:
    print "setting thk_initial = 0 ..."
    thk_init = 0.0 * thk_obs

nc.close()

# convert to PETSc-type vecs
x1vec = x1.view(pbio.Vec)
y1vec = y1.view(pbio.Vec)
topgvec = topg.flatten().view(pbio.Vec)
cmbvec = cmb.flatten().view(pbio.Vec)
thk_obsvec = thk_obs.flatten().view(pbio.Vec)
thk_initvec = thk_init.flatten().view(pbio.Vec)

# open petsc binary file
io = pbio.PetscBinaryIO()

# write fields **in a particular order**  (names do not matter)
print "writing vars x1,y1,topg,cmb,thk_obs,thk_init into %s ..." % args.outname
io.writeBinaryFile(args.outname, [x1vec,y1vec,topgvec,cmbvec,thk_obsvec,thk_initvec,])

