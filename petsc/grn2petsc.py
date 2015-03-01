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
args = parser.parse_args()

try:
    nc = NC(args.inname, 'r')
except:
    print "ERROR: can't read from file %s ..." % args.inname
    sys.exit(11)

x = nc.variables['x1'][:]
print "length x = %d" % (np.shape(x)[0])
y = nc.variables['y1'][:]
print "length y = %d" % (np.shape(y)[0])

# load data
topg = np.squeeze(nc.variables['topg'][:]).flatten()
print "length topg (flattened)                  = %d" % (np.shape(topg)[0])
cmb = np.squeeze(nc.variables['climatic_mass_balance'][:]).flatten()
print "length climatic_mass_balance (flattened) = %d" % (np.shape(topg)[0])

xvec = x.view(pbio.Vec)
yvec = y.view(pbio.Vec)
topgvec = topg.view(pbio.Vec)
cmbvec = cmb.view(pbio.Vec)

# open petsc binary file
io = pbio.PetscBinaryIO()

# write fields in a particular order; names won't actually matter
io.writeBinaryFile(args.outname, [xvec,yvec,topgvec,cmbvec,])

