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
parser.add_argument('--ratios', action="store_true",
                    help='generate cseps and ttime ratios between rsls and ssls; requires equal numbers')
args = parser.parse_args()

def readfile(filename):
    try:
        vfile = open(filename, 'r')
    except IOError:
        print 'ERROR: cannot open file %s' % filename
        sys.exit(1)
    x = vfile.readlines()
    vfile.close()
    if np.mod(len(x),6) != 0:
        print 'ERROR: file does not contain 6X lines'
        sys.exit(2)
    dx = []
    cseps = []
    maxD = []
    ttime = []
    datatype = []
    vitype = []
    for k in range(len(x)/6):
        dx.append(      float(x[6*k].rstrip()  ))
        cseps.append(   float(x[6*k+1].rstrip()))
        maxD.append(    float(x[6*k+2].rstrip()))
        ttime.append(   float(x[6*k+3].rstrip()))
        datatype.append(x[6*k+4].rstrip())  # sea or mcb
        vitype.append(  x[6*k+5].rstrip())  # rsls or ssls
    return dx, cseps, maxD, ttime, datatype, vitype

dx, cseps, maxD, ttime, datatype, vitype = readfile(args.inname)

# indices where grid is not coarsest:
nci = [i for i,h in enumerate(dx) if h < 7000.0]
# indices where rs or ss resp.
rsi = [i for i,x in enumerate(vitype) if ((x=='rsls') & (i in nci))]
ssi = [i for i,x in enumerate(vitype) if ((x=='ssls') & (i in nci))]
if args.ratios:
    if len(ssi) != len(rsi):
        print 'generating ratios requires equal numbers of rsls and ssls runs'
        sys.exit(99)
    print 'dx values:'
    print np.array([dx[i] for i in rsi])
    print 'continuation eps ratios ssls/rsls (big means ssls worse convergence):'
    erat = np.array([cseps[i] for i in ssi]) / np.array([cseps[i] for i in rsi])
    print erat
    print 'mean = %.3f' % np.mean(erat)
    print 'time ratios ssls/rsls (big means ssls bad):'
    trat = np.array([ttime[i] for i in ssi]) / np.array([ttime[i] for i in rsi])
    print trat
    print 'mean = %.3f' % np.mean(trat)
    sys.exit(0)

plt.figure(figsize=(7,5))
plt.hold(True)
for k in rsi+ssi:
    if vitype[k] == 'rsls':
        if datatype[k] == 'sea':
            mark = 'ok' # circle
            marksize = 10.0
        else:
            mark = 'sk'
            marksize = 9.0
    else:
        if datatype[k] == 'sea':
            mark = 'xk' # square
            marksize = 11.0
        else:
            mark = '+k'
            marksize = 13.0
    line, = plt.loglog(dx[k],cseps[k],mark,markersize=marksize,
                       markerfacecolor='w',markeredgewidth=1.5)
    if dx[k] > 4000.0:
        if datatype[k] == 'sea':
            DT = 'BM1'
        else:
            DT = 'MCB'
        if vitype[k] == 'rsls':
            VI = 'RS'
        else:
            VI = 'SS'
        line.set_label(DT + ' ' + VI)
dxmin, dxmax = 380.0, 6250.0
# lines at \eps_0 ... \eps_11:
klist = np.arange(12)
epslist = 10.0**(-klist/3.0)
for k in klist:
    plt.loglog([510.0,dxmax],[epslist[k],epslist[k]],'k--',lw=0.7)
    plt.text(420.0,0.9*epslist[k],r'$\epsilon_{%d}$' % k, fontsize=16.0)
plt.hold(False)
plt.axis([dxmin, dxmax, 0.0001/1.2, 2.0*1.0])
plt.xlabel(r'$\Delta x$', fontsize=16.0, labelpad=1)
plt.xticks([500, 1000, 2500, 5000], ('500m', '1km', '2.5km', '5km'), fontsize=12.0)
klist = np.arange(5)
plt.yticks(10.0**(-klist), ('1.0','0.1','0.01','0.001','0.0001'), fontsize=12.0)
ax = plt.axes()
ax.xaxis.grid(True)
ax.yaxis.grid(False)
plt.legend(fontsize=14.0, loc='upper right')
if args.o == None:
    plt.show()
else:
    plt.savefig(args.o,bbox_inches='tight')

