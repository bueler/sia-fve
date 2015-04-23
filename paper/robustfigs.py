#!/usr/bin/env python
# (C) 2015 Ed Bueler

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=
'''Generate robustness-of-Newton-convergence figures from ASCII files.''')
parser.add_argument('-o', metavar='FILENAME',
                    help='image filename (e.g. .png or .pdf)')
parser.add_argument('inname', metavar='FILE', type=str,
                    help='name of ASCII input file')
args = parser.parse_args()

def readfile(filename):
    try:
        vfile = open(filename, 'r')
    except IOError:
        print 'ERROR: cannot open file %s' % filename
        sys.exit(1)
    x = vfile.readlines()
    vfile.close()
    if np.mod(len(x),4) != 0:
        print 'ERROR: file does not contain 4X lines'
        sys.exit(2)
    res = []
    for k in range(len(x)/4):
        res.append((float(x[4*k].rstrip()), \
                    float(x[4*k+1].rstrip()), \
                    x[4*k+2].rstrip(), \
                    x[4*k+3].rstrip()))
    return res

res = readfile(args.inname)

#print res

plt.figure(figsize=(7,5))
plt.hold(True)
dx = []
eps = []
for k in range(len(res)):
    dx.append(res[k][0])
    eps.append(res[k][1])
    if res[k][2] == 'sea':
        mfc = 'w'
    else:
        mfc = 'k'
    if res[k][3] == 'rsls':
        mark = 'ok'
    else:
        mark = 'sk'
    line, = plt.loglog(dx[k],eps[k],mark,markersize=10.0,
                       markerfacecolor=mfc,markeredgewidth=2.0)
    if dx[k] > 7000.0:
        line.set_label(res[k][2]+' '+res[k][3])
# show \eps_0 ... \eps_11
dxmin, dxmax = 300.0, 12000.0
klist = np.arange(12)
epslist = 10.0**(-klist/3.0)
for k in klist:
    plt.loglog([370.0,dxmax],[epslist[k],epslist[k]],'k--',lw=0.7)
    plt.text(320.0,0.95*epslist[k],r'$\epsilon_{%d}$' % k)
plt.hold(False)
plt.axis([dxmin, dxmax, 0.0001/1.2, 2.0*1.0])
plt.xlabel(r'$\Delta x$',fontsize=16.0, labelpad=-3)
plt.xticks([500, 1000, 2500, 5000, 10000], ('500m', '1km', '2.5km', '5km', '10km'))
#plt.ylabel(r'$\epsilon$',fontsize=16.0)
klist = np.arange(5)
plt.yticks(10.0**(-klist), ('1','0.1','0.01','0.001','0.0001'))
#plt.grid(True)
ax = plt.axes()
ax.xaxis.grid(True)
ax.yaxis.grid(False)
plt.legend(fontsize=16.0)
if args.o == None:
    plt.show()
else:
    plt.savefig(args.o,bbox_inches='tight')

