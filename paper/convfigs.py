#!/usr/bin/env python
#
# (C) 2015 Ed Bueler

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Generate convergence figures from ASCII verification files.')
parser.add_argument('infile', metavar='FILENAME',
                    help='filename with verification results (e.g. dome.verif, bedstep.verif)')
parser.add_argument('-o', metavar='FILENAME',
                    help='image filename (e.g. .png or .pdf)')
parser.add_argument("--bedstep", action="store_true",
                    help="assume infile has results from bedrock step solution")
parser.add_argument("--bedsteptable", action="store_true",
                    help="generate relative volume table w Jarosch et al bedrock step results")
parser.add_argument("--showfile", action="store_true",
                    help="print contents of input file")
args = parser.parse_args()

def readveriffile(filename):
    try:
        vfile = open(filename, 'r')
    except IOError:
        print 'cannot open file ', filename
        sys.exit(1)
    linecount = 0;
    for line in vfile:
        if linecount == 0:
            N = (int)(line)
            A = np.zeros((6*N,))
        else:
            A[linecount-1] = (float)(line)
        linecount = linecount + 1
    A = np.reshape(A,(6,N)).T
    M = (linecount-1) / N
    A = A[:,:M]
    return N, A

N, A = readveriffile(args.infile)

if args.showfile:
    print N
    print A
    sys.exit(0)

dx = A[:,0]

# if want table, just do that
if args.bedsteptable:
    relvol = (A[:,4] - A[:,5]) / A[:,5]
    JSA = { 1000.0 : (-7.588, 116.912),
             500.0 : (-5.075, 132.095),
             250.0 : (-3.401, 139.384),
             125.0 : (-2.579, 142.997) }
    for j in range(len(dx)):
        if dx[j] in [1000.0, 500.0, 250.0, 125.0]:
            print "%4d & %6.3f & %6.3f & %6.3f \\\\" \
                % (dx[j],100.0*relvol[j],JSA[dx[j]][0],JSA[dx[j]][1])
    sys.exit(0)

# proceed with generating verification figure
maxerr = A[:,2]
averr  = A[:,3]
cav    = np.polyfit(np.log(dx),np.log(averr),1)
pav    = np.poly1d(cav)

plt.figure(figsize=(7,5))
plt.loglog(dx,maxerr,'ok',label='maximum error',markersize=10.0,
           markerfacecolor='w',markeredgewidth=1.0)
plt.hold(True)
plt.loglog(dx,averr,'*k',label='average error',markersize=10.0)
plt.loglog(dx,np.exp(pav(np.log(dx))),'--k',label=r'$O(\Delta x^{%.3f})$' % cav[0])
plt.hold(False)
plt.grid(True)
plt.ylabel('thickness error  (m)',fontsize=12.0)
if args.bedstep:
    plt.axis([100.0, 3000.0, 5.0, 200.0])
    plt.xticks(np.array([125.0, 250.0, 500.0, 1000.0, 2000.0]), ('125','250','500','1000','2000'))
    plt.xlabel(r'$\Delta x$  (m)',fontsize=12.0)
    plt.legend(fontsize=12.0,loc='upper right')
else:
    plt.axis([1.0e3, 150.0e3, 0.5, 2000.0])
    plt.xticks(np.array([2.0e3, 5.0e3, 10.0e3, 25.0e3, 50.0e3, 100.0e3]), ('2','5','10','25','50','100'))
    plt.xlabel(r'$\Delta x$  (km)',fontsize=12.0)
    plt.legend(fontsize=12.0,loc='lower right')

if args.o == None:
    plt.show()
else:
    plt.savefig(args.o,bbox_inches='tight')

