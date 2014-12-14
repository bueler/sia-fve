#!/usr/bin/env python
#
# (C) 2014 Ed Bueler
#
# Generate .png figures from ASCII VTK time-step files written by layer.c.
#
# Example of usage:
#   $ (cd ../petsc/ && make layer)
#   $ mkdir foo/
# Generate 100 time-steps:
#   $ ../petsc/layer -lay_genfigs foo/ -da_refine 3 -lay_dt 0.01 -lay_steps 100 -lay_adscheme 1
# Only turn every fourth into a .png:
#   $ python genfigs.py foo/ --stride 4
#   $ eog foo/*.png

import numpy
import argparse
import sys
from pylab import *

commandline = " ".join(sys.argv[:])

parser = argparse.ArgumentParser(description='FIXME')
# positional
parser.add_argument('prefix', metavar='PREFIX',
                    help='root of file names f-0-0.txt and u-$N-$T.txt',
                    default='')
# parameter
parser.add_argument('--stride', action='store', metavar='N', help='use slices 1:N:end',
                    default=1)
parser.add_argument('--uscale', action='store', metavar='N', help='scale thickness by this',
                    default=12.0)
args = parser.parse_args()
stride = int(args.stride)
uscale = float(args.uscale)

def readvecvtk(vfile,vfilename):
    headersread = 0
    count = 0
    for line in vfile:
        if line: # only act if line content remains
            # read headers
            if headersread == 0:
                vheader = line
                tmpstr = vheader.split()[0]
                N = (int)(vheader.split()[1])
                print 'reading N = %d values from %s ...' % (N,vfilename)
                v = numpy.zeros(N)
                headersread += 1
                continue
            elif (headersread == 1) or (headersread == 2):
                vheader = line
                headersread += 1
                continue
            elif (headersread >= 3) and (count >= N):
                break # nothing more to read
            # read content
            vline = line.split()
            v[count] = float(vline[0])
            count += 1
    return v, N

L = 10.0
def genfig(N, x, u, f, b):
    s = b + uscale * u
    plot(x, b, '--k', lw=3.0)
    hold(True)
    plot(x, s, 'k', lw=3.0)
    axis([0.0,L,-0.5,4.5])
    axis('off')
    astart = N/40
    astep  = N/20
    ascale = f.max() - f.min()
    for j in range(20):
        jj = astart + astep * j
        xarr = x[jj]
        magarr = f[jj] / ascale
        arrow(xarr,3.75,0.0,-magarr,lw=1.0,head_width=0.1,color='b')
    hold(False)

figdebug = False
def figsave(name):
    if figdebug:
        show()  # debug
    else:
        savefig(name,bbox_inches='tight')


for vec in ['x','f','b']:
    fname = args.prefix + vec + '.txt'
    try:
        ffile = open(fname, 'r')
    except IOError:
        print 'cannot open file ', fname
        sys.exit(1)
    if vec == 'x':
      x, _ = readvecvtk(ffile,fname)
    elif vec == 'f':
      f, _ = readvecvtk(ffile,fname)
    elif vec == 'b':
      b, _ = readvecvtk(ffile,fname)
    else:
      print 'how get here?'
      sys.exit(99)
    ffile.close()

huge = 1000  # limit number of input files, thus frames
pngk = 1
for txtk in range(1,huge,stride):
    uname = args.prefix + 'u-%d.txt' % txtk
    try:
        ufile = open(uname, 'r')
    except IOError:
        print 'cannot open file %s ... so done reading' % uname
        break
    u, Nu = readvecvtk(ufile,uname)
    ufile.close()
    ufig = figure(figsize=(10,4))
    genfig(Nu, x, u, f, b)
    figname = args.prefix + 'u-%d.png' % pngk
    figsave(figname)
    print '  figure file %s generated' % figname
    close(ufig)
    pngk = pngk + 1

