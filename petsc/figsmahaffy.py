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
# positional
#parser.add_argument('prefix', metavar='PREFIX',
#                    help='root of file names f-0-0.txt and u-$N-$T.txt',
#                    default='')
args = parser.parse_args()

def readvecvtk(vfile,vfilename):
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
    return v, NN

for vec in ['x', 'H', 'm', 'Herror']:
    fname = vec + '.txt'
    try:
        ffile = open(fname, 'r')
    except IOError:
        print 'cannot open file ', fname
        sys.exit(1)
    if vec == 'x':
      x, N = readvecvtk(ffile,fname)
    elif vec == 'H':
      H, NH = readvecvtk(ffile,fname)
    elif vec == 'm':
      m, Nm = readvecvtk(ffile,fname)
    elif vec == 'Herror':
      Herror, NHerror = readvecvtk(ffile,fname)
    else:
      print 'how get here?'
      sys.exit(99)
    ffile.close()

if (NH != Nm) | (NH != NHerror):
    print 'ERROR: different numbers of values in files'
    sys.exit(98)
if (NH != N*N):
    print 'ERROR: number of values in files is not a perfect square'
    sys.exit(98)

H = np.reshape(H,(N,N))
m = np.reshape(m,(N,N))
Herror = np.reshape(Herror,(N,N))

fsize = (12,9)

figdebug = False
def figsave(name):
    if figdebug:
        print '  showing %s ... close window to proceed ...' % name
        plt.show()  # debug
    else:
        plt.savefig(name,bbox_inches='tight')
        print '  figure file %s generated ...' % name

x = x/1000.0

plt.figure(figsize=fsize)
plt.pcolor(x,x,H)
plt.axis('tight')
plt.colorbar()
plt.title('thickness solution H  (m)')
figsave('H.png')

plt.figure(figsize=fsize)
plt.pcolor(x,x,m * 31556926.0)
plt.axis('tight')
plt.colorbar()
plt.title('surface mass balance m  (m/a)')
figsave('m.png')

plt.figure(figsize=fsize)
plt.pcolor(x,x,Herror)
plt.axis('tight')
plt.colorbar()
plt.title('thickness error H - H_exact  (m)')
figsave('Herror.png')

