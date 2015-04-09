#!/usr/bin/env python
# (C) 2015 Ed Bueler

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=
'''Generate convergence figures from ASCII verification files.
 Supply at least one input file name using options -bedstep or -dome.''')

parser.add_argument('-o', metavar='FILENAME', default='unnamed.pdf',
                    help='image filename (e.g. .png or .pdf)')
parser.add_argument('-bedstep', metavar='FILENAME', default='',
                    help="text file with results from bedrock step verification")
parser.add_argument('--bedsteptable', action="store_true",
                    help="generate relative volume table w Jarosch et al bedrock step results")
parser.add_argument('-tablefiles', default='', metavar='A,B', type=str,
                    help='comma-delimited list of files from which to get bedrock step results FOR TABLE')
parser.add_argument('-dome', metavar='FILENAME', default='',
                    help="text file with results from bedrock step verification")
parser.add_argument("--showdata", action="store_true",
                    help="print contents of input files")
parser.add_argument('-true', metavar='FILENAME', default='',
                    help="add results from dome verification using true Mahaffy; plotted on top of M* results")

args = parser.parse_args()

def readveriffile(filename):
    try:
        vfile = open(filename, 'r')
    except IOError:
        print 'cannot open file %s' % filename
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
    if args.showdata:
        print A.shape[0]
        print A
    return A

# if want table, just do that and exit
if args.bedsteptable:
    tabfiles = args.tablefiles.split(',')
    if len(tabfiles) != 3:
        print "ERROR: --bedsteptable is set but there is no LENGTH 3 -tablefiles filename list"
        sys.exit(2)
    epss = []
    relvols = []
    for j in range(len(tabfiles)):
        A = readveriffile(tabfiles[j])
        epss.append( A[:,1] )
        relvols.append( 100.0 * (A[:,4] - A[:,5]) / A[:,5] )
    dx = A[:,0]
    JSA = { 1000.0 : (-7.588, 116.912),
             500.0 : (-5.075, 132.095),
             250.0 : (-3.401, 139.384),
             125.0 : (-2.579, 142.997) }
    for j in range(len(dx)):
        if dx[j] in [1000.0, 500.0, 250.0, 125.0]:
            print "%4d & " % dx[j],
            for k in range(3):
                if epss[k][j] > 1.0e-3:
                    print "NC & ",   # not converged case
                else:
                    print "%6.3f & " % relvols[k][j],
            print "%6.3f & %6.3f \\\\" % (JSA[dx[j]][0],JSA[dx[j]][1])
    sys.exit(0)

infilelist = [args.bedstep, args.dome]
if (sum(map(len,infilelist)) == 0) or (np.prod(map(len,infilelist)) > 0):
    print "ERROR: exactly one of options -dome or -bedstep must be given with a file name"
    sys.exit(1)

def extract(A):
    return A[:,0], A[:,2], A[:,3], A[:,1]

# generate verification plots as image
for infile in infilelist:
    if len(infile)==0:
        continue
    dx, maxerr, averr, lasteps = extract(readveriffile(infile))
    notconv = (lasteps > 0.0)
    plt.figure(figsize=(7,5))
    plt.hold(True)
    plt.loglog(dx[notconv],maxerr[notconv],'ok',markersize=16.0,
               markerfacecolor='w',alpha=0.3,markeredgewidth=0.5)
    plt.loglog(dx,maxerr,'*k',label=r'max. error M$^\bigstar$',markersize=10.0,
               markerfacecolor='w',markeredgewidth=1.0)
    plt.loglog(dx[notconv],averr[notconv],'ok',markersize=16.0,
               markerfacecolor='w',alpha=0.3,markeredgewidth=0.5)
    plt.loglog(dx,averr,'*k',label=r'av. error M$^\bigstar$',markersize=8.0)
    if args.dome:
        cav = np.polyfit(np.log(dx[:5]),np.log(averr[:5]),1)
        pav = np.poly1d(cav)
        plt.loglog(dx[:5],np.exp(pav(np.log(dx[:5]))),':k')
        ratelabel=r'$O(\Delta x^{%.2f})$' % cav[0]
        plt.text(dx[2],0.5*averr[2],ratelabel,fontsize=16.0)

    if len(args.true)>0:
        dx, maxerr, averr, lasteps = extract(readveriffile(args.true))
        notconv = (lasteps > 0.0)
        plt.loglog(dx[notconv],maxerr[notconv],'ok',markersize=16.0,
                   markerfacecolor='w',alpha=0.3,markeredgewidth=0.5)
        plt.loglog(dx,maxerr,'sk',label='max. error Mahaffy',markersize=7.0,
                   markerfacecolor='w',markeredgewidth=1.0)
        plt.loglog(dx[notconv],averr[notconv],'ok',markersize=16.0,
                   markerfacecolor='w',alpha=0.3,markeredgewidth=0.5)
        plt.loglog(dx,averr,'sk',label='av. error Mahaffy',markersize=7.0)
        cav = np.polyfit(np.log(dx[:5]),np.log(averr[:5]),1)
        pav = np.poly1d(cav)
        plt.loglog(dx[:5],np.exp(pav(np.log(dx[:5]))),':k')
        ratelabel=r'$O(\Delta x^{%.2f})$' % cav[0]
        plt.text(dx[3],3.0*averr[3],ratelabel,fontsize=16.0)

    plt.hold(False)
    plt.grid(True)
    plt.ylabel('thickness error  (m)',fontsize=12.0)
    if args.bedstep:
        plt.axis([40.0, 3000.0, 5.0, 200.0])
        plt.xticks(np.array([62.5, 125.0, 250.0, 500.0, 1000.0, 2000.0]), ('62.5', '125','250','500','1000','2000'))
        plt.xlabel(r'$\Delta x$  (m)',fontsize=12.0)
        plt.legend(fontsize=12.0,loc='upper right')
    else:
        plt.axis([1.0e3, 150.0e3, 0.075, 2000.0])
        plt.xticks(np.array([2.0e3, 5.0e3, 10.0e3, 25.0e3, 50.0e3, 100.0e3]), ('2','5','10','25','50','100'))
        plt.xlabel(r'$\Delta x$  (km)',fontsize=12.0)
        plt.legend(fontsize=12.0,loc='lower right')

    if args.o == None:
        plt.show()
    else:
        plt.savefig(args.o,bbox_inches='tight')

