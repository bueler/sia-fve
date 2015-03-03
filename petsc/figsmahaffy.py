#!/usr/bin/env python
#
# (C) 2015 Ed Bueler
#
# Generate .png figures from PETSc binary files written by mahaffy.c.
#
# Example of usage:
#   $ make mahaffy
#   $ ./mahaffy -mah_dump
#   $ python figsmahaffy.py
#   $ eog *.png

import argparse
import sys

import numpy as np
import matplotlib.pyplot as plt

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

parser = argparse.ArgumentParser(description='Generate .png figures from PETSc binary files written by mahaffy.c.')
args = parser.parse_args()

def readvec(vfilename):
    try:
        FIXME
    except IOError:
        print 'cannot open file %s ... skipping ...' % vfilename
        return 0, -1

    return v, NN

for vec in ['x', 'y', 'H', 'b', 'm', 'Herror']:
    fname = vec + '.txt'
    if vec == 'x':
      x, Nx = readvec(fname)
    elif vec == 'y':
      y, Ny = readvec(fname)
    elif vec == 'H':
      H, NH = readvec(fname)
    elif vec == 'b':
      b, Nb = readvec(fname)
    elif vec == 'm':
      m, Nm = readvec(fname)
    elif vec == 'Herror':
      Herror, NHerror = readvec(fname)
    else:
      print 'how did I get here?'
      sys.exit(95)

if (NH != Nx*Ny):
    print 'ERROR: number of values in H.dat does not match axes in x.dat, y.dat'
    sys.exit(96)
H = np.reshape(H,(Ny,Nx))

if (NH != Nb):
    print 'ERROR: different sizes of fields H and b in files'
    sys.exit(97)
b = np.reshape(b,(Ny,Nx))

if (NH != Nm):
    print 'ERROR: different sizes of fields H and m in files'
    sys.exit(97)
m = np.reshape(m,(Ny,Nx))

if NHerror > 0:
    if (NH != NHerror):
        print 'ERROR: different sizes of fields H and Herror'
        sys.exit(99)
    Herror = np.reshape(Herror,(Ny,Nx))

figdebug = False
def figsave(name):
    if figdebug:
        print '  showing %s ... close window to proceed ...' % name
        plt.show()  # debug
    else:
        plt.savefig(name,bbox_inches='tight')
        print '  figure file %s generated ...' % name

x = x/1000.0
y = y/1000.0

fsize = (12,9)

plt.figure()
#plt.figure(figsize=fsize)
plt.pcolormesh(x,y,H)
plt.axis('equal')
plt.colorbar()
plt.title('thickness solution H (m) with %.2f <= H <= %.2f' % (H.min(),H.max()))
figsave('H.png')

plt.figure()
plt.pcolormesh(x,y,b)
plt.axis('equal')
if (b.max() > b.min()):
    plt.colorbar()
plt.title('bed elevation b (m) with %.2f <= b <= %.2f' % (b.min(),b.max()))
figsave('b.png')

plt.figure()
s = np.maximum(0.0, H + b)
plt.pcolormesh(x,y,s)
plt.axis('equal')
plt.colorbar()
plt.title('surface elevation s (m) with %.2f <= s <= %.2f' % (s.min(),s.max()))
figsave('s.png')

plt.figure()
m = m * 31556926.0
plt.pcolormesh(x,y,m)
plt.axis('equal')
plt.colorbar()
plt.title('surface mass balance m (m/a) with %.2f <= m <= %.2f' % (m.min(),m.max()))
figsave('m.png')

if NHerror > 0:
    plt.figure()
    plt.pcolormesh(x,y,Herror)
    plt.axis('equal')
    plt.colorbar()
    plt.title('thickness error Herror = H - Hexact (m) with %.2f <= Herror <= %.2f' % (Herror.min(),Herror.max()))
    figsave('Herror.png')

