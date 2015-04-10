#!/usr/bin/env python
# (C) 2015 Ed Bueler

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Generate Newton convergence figure from an ASCII files.')

parser.add_argument('inname', metavar='INNAME', default='',
                    help='name of ascii file with newton conv data (see Makefile)')
parser.add_argument('-o', metavar='OUTNAME',
                    help='image filename (e.g. .pdf)')

args = parser.parse_args()

def readfile(filename):
    try:
        vfile = open(filename, 'r')
    except IOError:
        print 'cannot open file %s' % filename
        sys.exit(1)
    A = []
    col = np.array([])
    for line in vfile:
        words = line.strip().split(' ')
        if 'CONVERGED' in words[0]+words[1]:
            A.append((0,col))
            col = np.array([])
        elif 'DIVERGED' in words[0]+words[1]:
            A.append((1,col))
            col = np.array([])
        elif (words[1] == 'SNES') and (words[2] == 'Function'):
            col = np.append(col,(float)(words[4]))
    return A

A = readfile(args.inname)

eps = 0.1**(np.arange(13)/3.0)
eps[-1] = 0.0

plt.figure(figsize=(7,5))
plt.hold(True)
convcount = 0
for j in range(len(A)):
    ishift = 0.8 + 1.2*j + 2.0*np.arange(len(A[j][1]))/len(A[j][1])
    plt.semilogy(ishift,A[j][1],'ok-',markersize=6.0,markerfacecolor='k')
    topres = A[j][1][0]
    plt.semilogy([ishift[0],ishift[0]], [1.5*topres,7.0*topres], 'k', lw=0.5)
    if A[j][0] == 0:
        tagfs = 16.0
        epsstr = r'$\epsilon_{%d}$' % convcount
        convcount += 1
    else:
        tagfs = 12.0
        epsstr = 'DIV'
    plt.text(ishift[0]-0.2,10.0*topres,epsstr,fontsize=tagfs)
plt.hold(False)
plt.grid(True,axis='y')
plt.ylabel('residual norm',fontsize=14.0)
plt.axis([0.0, ishift.max()+0.5, 2.0e-12, 1.0e5])
plt.xticks([], ())
#plt.xlabel('i',fontsize=16.0)
#plt.legend(fontsize=12.0,loc='upper right')

if args.o == None:
    plt.show()
else:
    plt.savefig(args.o,bbox_inches='tight')

