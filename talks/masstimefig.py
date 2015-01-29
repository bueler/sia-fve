#!/usr/bin/env python
#
# (C) 2015 Ed Bueler
#
# Generate .png figures from ASCII time-series mass files with t_n,M_n,R_n, as
# written by layer.c.
#
# Example of usage:
#   $ (cd ../petsc/ && make layer)
# Generate time-series in ASCII files:
#   $ ../petsc/layer -lay_dt 0.02 -lay_steps 500 -lay_timedependentsource -lay_massfile foo1.txt
#   $ ../petsc/layer -lay_dt 0.01 -lay_steps 1000 -lay_timedependentsource -lay_massfile foo2.txt
# Make image file from time-series:
#   $ python masstimefig.py -o bar.pdf foo1.txt foo2.txt
# or simply
#   $ python masstimefig.py -o bar.pdf foo?.txt

import numpy
import argparse
import sys
from pylab import *

parser = argparse.ArgumentParser(description='Generate .png figures from time-series mass files written by layer.c.')
parser.add_argument('files', metavar='FILE', nargs='+',
                    help='an ASCII file with these numbers in its first three columns: t_n M_n R_n')
parser.add_argument('-o', metavar='OUTFILE',
                    help='name of output image file (e.g. .png or .pdf); uses show() if not given',
                    default='')
args = parser.parse_args()

def readtimeseries(filename):
    S = loadtxt(filename)
    return S[:,0], S[:,1], S[:,2]

for name in args.files:
    t, M, R = readtimeseries(name)
    semilogy(t,M,'k')
    hold(True)
    semilogy(t,R,label='dt=%.3f' % (t[1]-t[0]))

hold(False)
axis([0.0,t.max(),1.0e-10,1.0])
legend(loc='lower right',fontsize=14.0)
text(4.0,0.2,r'$M_n$',fontsize=18.0)
text(4.9,1.0e-5,r'$R_n$',fontsize=18.0)
xlabel('t',fontsize=12.0)
ylabel('volume',fontsize=12.0)

if len(args.o) > 0:
    savefig(args.o,bbox_inches='tight')
else:
    show()

