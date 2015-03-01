#!/usr/bin/env python
#
# (C) 2015 Ed Bueler
#
# Generate .png figures from ASCII VTK files written by mahaffy.c.
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

parser = argparse.ArgumentParser(description='Generate .png figures from ASCII VTK files written by mahaffy.c.')
args = parser.parse_args()

def readvecvtk(vfilename):
    try:
        vfile = open(vfilename, 'r')
    except IOError:
        print 'cannot open file %s ... skipping ...' % vfilename
        return 0, -1
    headersread = 0
    count = 0
    for line in vfile:
        if line: # only act if line content remains
            # read headers
            if headersread == 0:
                vheader = line
                tmpstr = vheader.split()[0]
                NN = (int)(vheader.split()[1])
                print '  reading N = %d values from %s ...' % (NN,vfilename)
                v = np.zeros(NN)
                headersread += 1
                continue
            elif (headersread == 1) or (headersread == 2):
                vheader = line
                headersread += 1
                continue
            elif (headersread >= 3) and (count >= NN):
                break # nothing more to read
            # read content
            vline = line.split()
            v[count] = float(vline[0])
            count += 1
    vfile.close()
    return v, NN

for vec in ['x', 'y', 'H', 'b', 'm', 'Herror']:
    fname = vec + '.txt'
    if vec == 'x':
      x, Nx = readvecvtk(fname)
    elif vec == 'y':
      y, Ny = readvecvtk(fname)
    elif vec == 'H':
      H, NH = readvecvtk(fname)
    elif vec == 'b':
      b, Nb = readvecvtk(fname)
    elif vec == 'm':
      m, Nm = readvecvtk(fname)
    elif vec == 'Herror':
      Herror, NHerror = readvecvtk(fname)
    else:
      print 'how get here?'
      sys.exit(95)

if (NH != Nx*Ny):
    print 'ERROR: number of values in H.txt does not match axes in x.txt, y.txt'
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
plt.title('thickness solution H  (m)')
figsave('H.png')

plt.figure()
plt.pcolormesh(x,y,b)
plt.axis('equal')
if (b.max() > b.min()):
    plt.colorbar()
plt.title('bed elevation b  (m)')
figsave('b.png')

if (b.max() > b.min()):
    plt.figure()
    usurf = np.maximum(0.0, H + b)
    plt.pcolormesh(x,y,usurf)
    plt.axis('equal')
    plt.colorbar()
    plt.title('surface elevation usurf  (m)')
    figsave('usurf.png')
else:
    print 'not generating usurf.png because bed is flat'

plt.figure()
plt.pcolormesh(x,y,m * 31556926.0)
plt.axis('equal')
plt.colorbar()
plt.title('surface mass balance m  (m/a)')
figsave('m.png')

if NHerror > 0:
    plt.figure()
    plt.pcolormesh(x,y,Herror)
    plt.axis('equal')
    plt.colorbar()
    plt.title('thickness error H - H_exact  (m)')
    figsave('Herror.png')

