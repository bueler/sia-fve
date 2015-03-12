#!/usr/bin/env python
#
# (C) 2015 Ed Bueler

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Generate convergence figures from ASCII verification files dome.verif, bedstep.verif.')
parser.add_argument('prefix', metavar='PREFIX/',
                    help='prefix for verification results files dome.verif, bedstep.verif')
parser.add_argument('-o', metavar='OUTFILE',
                    help='name of output file (e.g. .png or .pdf); uses show() if not given',
                    default='')
args = parser.parse_args()

def readveriffile(filename):
    vname = args.prefix + filename
    try:
        vfile = open(vname, 'r')
    except IOError:
        print 'cannot open file ', vname
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

#for name in ['dome.verif','bedstep.verif']:
#    N, A = readveriffile(name)
#    print N
#    print A

Ndome, Adome = readveriffile('dome.verif')
dx     = Adome[:,0]
maxerr = Adome[:,2]
averr  = Adome[:,3]
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
plt.axis([1.5e3, 200.0e3, 0.5, 2000.0])
plt.xticks(np.array([2.0e3, 5.0e3, 1.0e4, 2.0e4, 5.0e4, 1.0e5]), ('2','5','10','20','50','100'))
plt.legend(fontsize=12.0,loc='lower right')
plt.ylabel('thickness error  (m)',fontsize=12.0)
plt.xlabel(r'$\Delta x$  (km)',fontsize=12.0)

# FIXME add bedstep fig

if len(args.o) > 0:
    plt.savefig(args.o,bbox_inches='tight')
else:
    plt.show()

