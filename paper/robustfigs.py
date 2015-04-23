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
    if np.mod(len(x),3) != 0:
        print 'ERROR: file does not contain 3X lines'
        sys.exit(2)
    res = []
    for k in range(len(x)/3):
        res.append((float(x[3*k].rstrip()), \
                    float(x[3*k+1].rstrip()), \
                    x[3*k+2].rstrip()))
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
    if res[k][2] == 'rsls':
        mark = 'ok'
    else:
        mark = 'sk'
    line, = plt.loglog(dx[k],eps[k],mark,markersize=10.0,
                       markerfacecolor='w',markeredgewidth=2.0)
    if k < 2:
        line.set_label(res[k][2])
plt.hold(False)
plt.axis([min(dx)/1.2, 1.2*max(dx), 0.5*min(eps), 1.2*max(eps)])
plt.grid(True)
plt.xlabel(r'$\Delta x$',fontsize=16.0, labelpad=-10)
plt.ylabel(r'$\epsilon$',fontsize=16.0)
plt.legend(fontsize=16.0, loc='lower right')
if args.o == None:
    plt.show()
else:
    plt.savefig(args.o,bbox_inches='tight')

