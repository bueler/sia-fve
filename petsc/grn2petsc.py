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

vec = np.array([1., 2., 3.]).view(PetscBinaryIO.Vec)
io = pbio.PetscBinaryIO()
io.writeBinaryFile('file.dat', [vec,])

parser = argparse.ArgumentParser(description='Generate PETSc binary format file from Greenland NetCDF file.')
# positional
parser.add_argument('inname', metavar='INNAME',
                    help='name of NetCDF input file with topg variable (e.g. pism_Greenland_5km_v1.1.nc)',
                    default='')
parser.add_argument('outname', metavar='OUTNAME',
                    help='name of output PETSc binary file (e.g. grn.dat)',
                    default='')
#parser.add_argument('--uscale', action='store', metavar='N', help='scale thickness by this',default=12.0)
#parser.add_argument('--taxis', action='store', metavar='N', type=int, help='if set, add time axis along bottom, marked at time-step, with N total timesteps',default=0)
args = parser.parse_args()


try:
    nc = NC(args.inname, 'r')
except:
    print "ERROR: can't read from file %s ..." % args.inname
    sys.exit(11)

x = nc.variables['x'][:]
y = nc.variables['y'][:]
width  = x.max() - x.min()
height = y.max() - y.min()
    

# load data
topg = np.squeeze(nc.variables['topg'][:])

topgvec = numpy.array([1., 2., 3.]).view(pbio.Vec)

io = pbio.PetscBinaryIO()
io.writeBinaryFile(args.outname, [topgvec,])


