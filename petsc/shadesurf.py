#!/usr/bin/env python
# (C) 2015 Ed Bueler

# example from http://matplotlib.org/1.4.2/examples/pylab_examples/shading_example.html

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource

import argparse

parser = argparse.ArgumentParser(description='Generate shaded figure of surface elevation from PETSc binary file written by mahaffy.c.')
parser.add_argument('-i', default='unnamed.dat', metavar='NAME', type=str,
                    help='name of special format PETSc binary file from which to read; default=%(default)s)')
args = parser.parse_args()

try:
    from petsc2numpy import readvecs
except:
    print "importing petsc2numpy failed"
    sys.exit(2)

print "reading from file %s ..." % args.i
d = readvecs(args.i,num=2,metaname='dimension variable',names=['x','y'])
x,y = d[0],d[1]
v = readvecs(args.i,num=7,skip=2,shape=(len(y),len(x)),names=['b','m','Hexact','H','residual'])
#b,_,_,H,_ = v[0],v[1],v[2],v[3],v[4]
b = v[0]
H = v[3]

print 'read b,H from file %s; shapes (%d,%d)' % (args.i,len(y),len(x))

def gets(H,b):
    Hdraft = -(910.0/1028.0) * H
    icebase = np.maximum(b,Hdraft)
    s = H + icebase
    return s

s = gets(H,b)
print 'surface computed ...'

# create light source object.
#ls = LightSource(azdeg=0,altdeg=65)
ls = LightSource(azdeg=0,altdeg=50)

# shade data, creating an rgb array.
rgb = ls.shade(np.flipud(s),plt.cm.gray)

print 'light source created ...'

plt.figure(figsize=(5,12))
plt.imshow(rgb)
print 'figure generated ...'

#plt.title('imshow with shading')
plt.xticks([])
plt.yticks([])
plt.show()

