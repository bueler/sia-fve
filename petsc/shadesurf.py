#!/usr/bin/env python
# (C) 2015 Ed Bueler

# example from http://matplotlib.org/1.4.2/examples/pylab_examples/shading_example.html

import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import argparse

parser = argparse.ArgumentParser(description='Generate shaded figure of surface elevation from PETSc binary file written by mahaffy.c.')
parser.add_argument('-o', metavar='FILENAME',
                    help='image filename (e.g. .png or .pdf)')
parser.add_argument('inname', metavar='FILE.dat', type=str,
                    help='name of .dat input file')
parser.add_argument('-dpi', default=200, type=int,
                    help='save figure at this resolution')
args = parser.parse_args()

try:
    from petsc2numpy import readvecs
except:
    print "importing petsc2numpy failed"
    sys.exit(2)

print "reading from file %s ..." % args.inname
d = readvecs(args.inname,num=2,metaname='dimension variable',names=['x','y'])
x,y = d[0],d[1]
v = readvecs(args.inname,num=7,skip=2,shape=(len(y),len(x)),names=['b','m','Hexact','H','residual'])
#b,_,_,H,_ = v[0],v[1],v[2],v[3],v[4]
b = v[0]
H = v[3]

print 'read b,H from file %s; shapes (%d,%d)' % (args.inname,len(y),len(x))

def gets(H,b):
    Hdraft = -(910.0/1028.0) * H
    icebase = np.maximum(b,Hdraft)
    s = H + icebase
    return s

s = gets(H,b)
print 'surface computed ...'

#alternative: see shading.py
#from shading import set_shade, hillshade
#rgb = set_shade(np.flipud(s),cmap=plt.cm.gist_earth)
#rgb = set_shade(np.flipud(s),cmap=plt.cm.gray,azdeg=300.0,altdeg=30.0)

# create light source object.
#ls = LightSource(azdeg=0,altdeg=65)
ls = LightSource(azdeg=300,altdeg=50)

# shade data, creating an rgb array.
rgb = ls.shade(np.flipud(s),plt.cm.gray)

print 'light source created ...'

plt.figure(figsize=(5,12))
plt.imshow(rgb, extent=(0,len(x),0,len(y)), interpolation='nearest')
#plt.imshow(rgb, interpolation='nearest')
plt.xticks(visible=False)
plt.yticks(visible=False)

ax = plt.axes()

axins = zoomed_inset_axes(ax, 5, loc=1) # zoom = 3
axins.imshow(rgb, extent=(0,len(x),0,len(y)), interpolation='nearest')
#axins.imshow(rgb, interpolation='nearest')

#x1, x2, y1, y2 = -1.5, -0.9, -2.5, -1.9
aspect = 1.5
ox = 0.20 * len(x)
oy = 0.30 * len(y)
width = 0.12 * len(x)
axins.set_xlim(ox, ox+width)
axins.set_ylim(oy, oy+aspect*width)

plt.xticks(visible=False)
plt.yticks(visible=False)

mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

print 'figure generated ...'

#plt.title('imshow with shading')
#plt.xticks([])
#plt.yticks([])

if args.o == None:
    plt.show()
else:
    plt.savefig(args.o,dpi=args.dpi,bbox_inches='tight')


