#!/usr/bin/env python
# (C) 2015 Ed Bueler
#
# Bilinearly-remap fields from SeaRISE file to grid used in mcbRES.nc (derived
# from MCdataset-2014-11-19.nc by running getmcb.sh):
# Specifically, remap:
#    in SeaRISE projection:                        in MCB projection:
#    climatic_mass_balance       --> remap -->     cmb
#    topg<0 (bathymetry)         --> merge -->     topg

import argparse
import sys
import numpy as np

try:
    from netCDF4 import Dataset as NC
except:
    print "netCDF4 (=netcdf4-python) is not installed? ... see http://unidata.github.io/netcdf4-python/"
    sys.exit(1)

try:
    from pyproj import Proj, transform
except:
    print "pyproj is not installed? ... see http://jswhit.github.io/pyproj/"
    sys.exit(2)

def get_dims(nc):
    # possible dimension names
    xdims = ['x', 'x1']
    ydims = ['y', 'y1']
    zdims = ['z', 'z1']
    tdims = ['t', 'time']
    xdim,ydim,zdim,tdim = None,None,None,None
    # assign dimensions
    for dim in xdims:
        if dim in list(nc.dimensions.keys()):
            xdim = dim
    for dim in ydims:
        if dim in list(nc.dimensions.keys()):
            ydim = dim
    for dim in zdims:
        if dim in list(nc.dimensions.keys()):
            zdim = dim
    for dim in tdims:
        if dim in list(nc.dimensions.keys()):
            tdim = dim
    return xdim, ydim, zdim, tdim

parser = argparse.ArgumentParser(description='Bilinearly-remap fields from SeaRISE file to grid used in MCdataset-2014-11-19.nc.')
parser.add_argument('inname', metavar='INNAME',
                    help='input NetCDF file (SeaRISE-type grid)',
                    default='')
parser.add_argument('targetname', metavar='TARGET',
                    help='NetCDF file with target grid (e.g. mcb4500m.nc)',
                    default='')
parser.add_argument('outname', metavar='OUTNAME',
                    help='output NetCDF file to create (e.g. fix4500m.nc)',
                    default='')
args = parser.parse_args()

try:
    ncsrc = NC(args.inname, 'r')
except:
    print "ERROR: can't open file %s for reading ..." % args.inname
    sys.exit(11)

try:
    nctarg = NC(args.targetname, 'r')
except:
    print "ERROR: can't open file %s for reading ..." % args.targetname
    sys.exit(12)

try:
    ncout = NC(args.outname, 'w', format='NETCDF3_CLASSIC')
except:
    print "ERROR: can't open file %s for writing ..." % args.outname
    sys.exit(13)

# read from source file
xdim, ydim, zdim, tdim = get_dims(ncsrc)
x1src = ncsrc.variables[xdim][:].astype(np.float64)  # is increasing
y1src = ncsrc.variables[ydim][:].astype(np.float64)  # is increasing
print "source grid from %s has dimensions (%s,%s)=(%d,%d)" % (args.inname,xdim,ydim,len(y1src),len(x1src))
print "    reading cmb and topg ..."
cmbsrc = np.squeeze(ncsrc.variables['cmb'][:].astype(np.float64))
topgsrc = np.squeeze(ncsrc.variables['topg'][:].astype(np.float64))

# NOTE  ncsrc.proj4 is type 'unicode', so convert to type 'str'
projsrc = Proj(ncsrc.proj4.encode('ascii','ignore'))

# read from target file
xdim, ydim, zdim, tdim = get_dims(nctarg)
xtarg = nctarg.variables[xdim][:].astype(np.float64)  # is increasing
ytarg = nctarg.variables[ydim][:].astype(np.float64)  # is *decreasing*
print "target grid from %s has dimensions (%s,%s)=(%d,%d)" % (args.targetname,xdim,ydim,len(ytarg),len(xtarg))
print "    reading thk and topg_nobathy ..."
thktarg = np.squeeze(nctarg.variables['thk'][:].astype(np.float64))
topgnobathtarg = np.squeeze(nctarg.variables['topg_nobathy'][:].astype(np.float64))

# NOTE  nctarg.proj4 = +init=espg:3413 is misspelled, so replace it with string literal
projtarg = Proj('+init=epsg:3413')

def insrc(x,y):
    return (x >= x1src.min()) & (x <= x1src.max()) & (y >= y1src.min()) & (y <= y1src.max())

def instr(flg):
    if flg:
        return "in source grid"
    else:
        return "outside"

def demoremap():
    # show where four corners and middle and more points of target will go in src region
    ix = [  0,   0, 333, 333, 166,  50,  50, 300, 300]
    iy = [  0, 598,   0, 598, 298,  50, 550,  50, 550]
    for k in range(len(ix)):
        x, y = transform(projtarg, projsrc, xtarg[ix[k]], ytarg[iy[k]])
        flg = insrc(x, y)
        print "(%11.2f,%11.2f) --> (%11.2f,%11.2f)  ... %s" % \
              (xtarg[ix[k]], ytarg[iy[k]], x, y, instr(flg))
        if flg:
            iix = np.nonzero(x1src <= x)[0].max()
            iiy = np.nonzero(y1src <= y)[0].max()
            print "    box is [%11.2f,%11.2f] x [%11.2f,%11.2f]" % \
                  (x1src[iix],x1src[iix+1],y1src[iiy],y1src[iiy+1])

#demoremap()

# fbox variable is traversed in counterclockwise order from lower left:
#     fbox[3]           fbox[2]
#   xx[0],yy[1]       xx[1],yy[1]
#        +-----------------+
#        |                 |
#        |                 |
#        +-----------------+
#   xx[0],yy[0]       xx[1],yy[0]
#     fbox[0]           fbox[1]
def bilin(x,y,xx,yy,fbox):
    if (xx[0]>=xx[1]) or (yy[0]>=yy[1]):
        raise ValueError("bilin has xx or yy which is not increasing")
    if (x<xx[0]) or (x>xx[1]) or (y<yy[0]) or (y>yy[1]):
        print "DEBUG: rectangle xx and yy are: ", xx, yy
        raise ValueError("bilin has (x,y)=(%f,%f) outside of rectangle" % (x,y))
    xi   = (x - xx[0]) / (xx[1] - xx[0])
    eta  = (y - yy[0]) / (yy[1] - yy[0])
    cxi  = 1.0 - xi
    ceta = 1.0 - eta
    return (fbox[0] * cxi + fbox[1] * xi) * ceta + (fbox[2] * xi + fbox[3] * cxi) * eta

def evalfbox(var,iix,iiy):
    return np.array([var[iiy][iix], var[iiy][iix+1], var[iiy+1][iix+1], var[iiy+1][iix]])

print "bilinearly remapping cmb and merging topg with bathymetry ..."
cmbtarg = np.zeros((len(ytarg),len(xtarg)))
topgtarg = np.zeros((len(ytarg),len(xtarg)))
cmb_fill = -10.0 / 31556926.0  # = 10.0 m year-1
topg_fill = -300.0  # m a.s.l.
for k in range(len(ytarg)):
    for j in range(len(xtarg)):
        x, y = transform(projtarg, projsrc, xtarg[j], ytarg[k])
        flg = insrc(x, y)
        if flg:
            iix = np.nonzero(x1src <= x)[0].max()
            iiy = np.nonzero(y1src <= y)[0].max()
            if (iix+1 < len(x1src)) and (iiy-1 >= 0):
                xx = x1src[iix:iix+2]
                yy = y1src[iiy:iiy+2]
            else:
                flg = False
        if flg:
            cmbtarg[k][j]  = bilin(x,y,xx,yy,evalfbox(cmbsrc,iix,iiy))
            if topgnobathtarg[k][j] < -9900.0:
                topgtarg[k][j] = bilin(x,y,xx,yy,evalfbox(topgsrc,iix,iiy))
            else:
                topgtarg[k][j] = topgnobathtarg[k][j]
        else:
            cmbtarg[k][j]  = cmb_fill
            topgtarg[k][j] = topg_fill
    if (k>0) and (np.mod(k,100)==0):
        print k,
    elif (k>0) and (np.mod(k,10)==0):
        print ".",
    sys.stdout.flush()
print " "

ncsrc.close()
nctarg.close()

def deftargvar(nc, name, units):
    var = nc.createVariable(name, 'f4', dimensions=("y1", "x1"))
    var.units = units
    return var

print "creating dimensions and variable x,y in %s ..." % args.outname
ncout.createDimension("x1", size=len(xtarg))
ncout.createDimension("y1", size=len(ytarg))
x_var = ncout.createVariable("x1", 'f4', dimensions=("x1",))
x_var.units = "m"
x_var.long_name = "easting"
x_var.standard_name = "projection_x_coordinate"
x_var[:] = xtarg
y_var = ncout.createVariable("y1", 'f4', dimensions=("y1",))
y_var.units = "m"
y_var.long_name = "northing"
y_var.standard_name = "projection_y_coordinate"
y_var[:] = ytarg[::-1]

mapping = ncout.createVariable("mapping", 'c')
mapping.grid_mapping_name = "polar_stereographic"
mapping.straight_vertical_longitude_from_pole = -45. ;
mapping.false_easting = 0. ;
mapping.false_northing = 0. ;
mapping.latitude_of_projection_origin = 90. ;
mapping.standard_parallel = 70. ;


print "copying thk variable into %s ..." % args.outname
thk_var = deftargvar(ncout, "thk", "m")
thk_var.standard_name = "land_ice_thickness"
thk_var.grid_mapping = "mapping"
thk_var[:] = np.flipud(thktarg)

print "putting new cmb variable in %s ..." % args.outname
cmb_var = deftargvar(ncout, "cmb", "m s-1")
cmb_var.grid_mapping = "mapping"
cmb_var.standard_name = "land_ice_surface_specific_mass_balance"
cmb_var[:] = np.flipud(cmbtarg)

print "putting new topg variable in %s ..." % args.outname
topg_var = deftargvar(ncout, "topg", "m")
topg_var.standard_name = "bedrock_altitude"
topg_var.grid_mapping = "mapping"
topg_var[:] = np.flipud(topgtarg)

ncout.Conventions = "CF-1.6"
ncout.proj4 = "+init=epsg:3413"
ncout.close()

