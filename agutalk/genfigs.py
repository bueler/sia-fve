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
#   $ python genfigs.py ../petsc/foo/ --stride 4
#   $ eog ../petsc/foo/*.png

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
# integer
parser.add_argument('--stride', action='store', metavar='N',
                    help='use slices 1:N:end',
                    default=1)
args = parser.parse_args()
stride = int(args.stride)

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
    vfile.close()
    return v, N

L = 10.0
def genfig(N, u, f):
    dx = L / float(N)
    x = linspace(0.0+dx/2.0,L-dx/2.0,N)
    b = 0.07*(x-3.0)**2 + 0.2*sin(2.0*x)  # see figs1d.py
    uscale = 12.0
    s = b + uscale * u
    plot(x, b, '--k', lw=3.0)
    hold(True)
    plot(x, s, 'k', lw=3.0)
    axis([0.0,L,-0.5,4.5])
    #axis('off')
    astart = N/40
    astep  = N/20
    ascale = f.max() - f.min()
    for j in range(20):
        jj = astart + astep * j
        xarr = x[jj]
        magarr = f[jj] / ascale
        arrow(xarr,3.0,0.0,magarr,lw=1.0,head_width=0.1,color='b')
    hold(False)
    return x, b

figdebug = False
def figsave(name):
    if figdebug:
        show()  # debug
    else:
        savefig(name,bbox_inches='tight')


fname = args.prefix + 'f.txt'
try:
    ffile = open(fname, 'r')
except IOError:
    print 'cannot open file ', fname
    sys.exit(1)
f, Nf = readvecvtk(ffile,fname)

huge = 1000  # limit number of frames
for ucount in range(1,huge,stride):
    uroot = args.prefix + 'u-%d' % ucount
    uname = uroot + '.txt'
    try:
        ufile = open(uname, 'r')
    except IOError:
        print 'cannot open file %s ... so done reading' % uname
        break
    u, Nu = readvecvtk(ufile,uname)
    ufig = figure(figsize=(10,4))
    _, _ = genfig(Nu, u, f)
    figname = '%s.png' % uroot
    figsave(figname)
    print '  figure file %s generated' % figname
    close(ufig)

