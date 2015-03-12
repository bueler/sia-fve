#!/usr/bin/env python
#
# (C) 2015 Ed Bueler

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Generate convergence figures from ASCII verification files.')
parser.add_argument('prefix', metavar='PREFIX/',
                    help='prefix for verification results files dome.verif, bedstep.verif')
parser.add_argument("--pdf", help="generate PDF files domeconv.pdf, bedstepconv.pdf, bedsteprelvol.pdf",
                    action="store_true")
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

def showcontents():
    for name in ['dome.verif','bedstep.verif']:
        N, A = readveriffile(name)
        print N
        print A

#showcontents()

Ndome, Adome = readveriffile('dome.verif')
dx     = Adome[:,0]
maxerr = Adome[:,2]
averr  = Adome[:,3]
cav    = np.polyfit(np.log(dx),np.log(averr),1)
pav    = np.poly1d(cav)

plt.figure(figsize=(7,5))
plt.loglog(dx/1000.0,maxerr,'ok',label='maximum error',markersize=10.0,
           markerfacecolor='w',markeredgewidth=1.0)
plt.hold(True)
plt.loglog(dx/1000.0,averr,'*k',label='average error',markersize=10.0)
plt.loglog(dx/1000.0,np.exp(pav(np.log(dx))),'--k',label=r'$O(\Delta x^{%.3f})$' % cav[0])
plt.hold(False)
plt.grid(True)
plt.axis([1.5, 200.0, 0.5, 2000.0])
plt.xticks(np.array([2.0, 5.0, 10.0, 20.0, 50.0, 100.0]), ('2','5','10','20','50','100'))
plt.legend(fontsize=12.0,loc='lower right')
plt.ylabel('thickness error  (m)',fontsize=12.0)
plt.xlabel(r'$\Delta x$  (km)',fontsize=12.0)
if args.pdf:
    plt.savefig('domeconv.pdf',bbox_inches='tight')

Nbs, Abs = readveriffile('bedstep.verif')
dx     = Abs[:,0]
maxerr = Abs[:,2]
averr  = Abs[:,3]
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
plt.axis([200.0, 7000.0, 5.0, 200.0])
plt.xticks(np.array([300.0, 1000.0, 5000.0]), ('300','1000','5000'))
plt.legend(fontsize=12.0,loc='lower right')
plt.ylabel('thickness error  (m)',fontsize=12.0)
plt.xlabel(r'$\Delta x$  (km)',fontsize=12.0)
if args.pdf:
    plt.savefig('bedstepconv.pdf',bbox_inches='tight')

relvol = abs(Abs[:,4] - Abs[:,5]) / Abs[:,5]
plt.figure(figsize=(7,5))
plt.loglog(dx,100.0*relvol,'*k',markersize=10.0)
plt.grid(True)
plt.axis([200.0, 7000.0, 0.05, 50.0])
plt.xticks(np.array([300.0, 1000.0, 5000.0]), ('300','1000','5000'))
plt.ylabel('relative volume error  (%)',fontsize=12.0)
plt.xlabel(r'$\Delta x$  (km)',fontsize=12.0)
if args.pdf:
    plt.savefig('bedsteprelvol.pdf',bbox_inches='tight')

if not args.pdf:
    plt.show()

