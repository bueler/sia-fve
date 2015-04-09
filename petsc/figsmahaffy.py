#!/usr/bin/env python
# (C) 2015 Ed Bueler
#
# Example of usage: see README.md and ../../mahaffypaper/Makefile.

import argparse
from exactforfigs import dome_exact, bedstep_exact

parser = argparse.ArgumentParser(description='Generate figures from PETSc binary files written by mahaffy.c.  Default is to generate both .pdf profile and .png map-plane figures.')
parser.add_argument("--profile", action="store_true",
                    help="generate profile figure only")
parser.add_argument("--half", action="store_true",
                    help="make half-width profile")
parser.add_argument("--observed", action="store_true",
                    help="label Hexact as observed in profile")
parser.add_argument("--exactbed", action="store_true",
                    help="plot exact solution from Jarosch et al (2013)")
parser.add_argument("--exactdome", action="store_true",
                    help="plot exact solution from Bueler (2003)")
parser.add_argument("--blowup", action="store_true",
                    help="plot detail inset; use only with --exactdome")
parser.add_argument("--map", action="store_true",
                    help="generate map-plane figures only")
parser.add_argument('-extra_H', default='', metavar='A,B', type=str,
                    help='comma-delimited list of files from which to ADDITIONALLY read thickness H, for combined profile')
parser.add_argument('-extra_H_labels', default='', metavar='"0","A","B"', type=str,
                    help='comma-delimited list of ALL labels for combined profile; length must be one more than -extra_H list')
args = parser.parse_args()

if (not args.profile) & (not args.map):
    args.profile = True
    args.map = True

extraH = (len(args.extra_H) > 0)
if extraH:
    if args.map:
        print "option -extra_H only allowed for profile plot"
        sys.exit(3)
    Hlist = args.extra_H.split(',')
    Hlabels = args.extra_H_labels.split(',')
    print Hlist
    print Hlabels
    if len(Hlist)+1 != len(Hlabels):
        print "length of labels list must be one more than -extra_H list"
        sys.exit(4)

import sys
import numpy as np
import matplotlib.pyplot as plt

try:
    import PetscBinaryIO as pbio
except:
    print "'import PetscBinaryIO' failed"
    print "need link to petsc/bin/petsc-pythonscripts/PetscBinaryIO.py?"
    sys.exit(2)

try:
    import petsc_conf
except:
    print "'import petsc_conf.py' failed"
    print "need link to petsc/bin/petsc-pythonscripts/petsc_conf.py?"
    sys.exit(2)

# read a vec from a petsc binary file
io = pbio.PetscBinaryIO()
def readvec(fname,shape=(0,0)):
    try:
        fh = open(fname)
    except:
        print "unable to open '%s' ... ending ..." % fname
        sys.exit(3)
    objecttype = io.readObjectType(fh)
    if objecttype == 'Vec':
        v = io.readVec(fh)
        fh.close()
        if sum(shape) > 0:
            v = np.reshape(v,shape)
        print "  read vec from '%s' with shape %s" % (fname, str(np.shape(v)))
        return v
    else:
        print "unexpected objectype '%s' ... ending ..." % objecttype
        sys.exit(4)

print "reading vars x,y,b,m,Hexact,H from *.dat ..."
x = readvec('x.dat')
y = readvec('y.dat')
b = readvec('b.dat',shape=(len(y),len(x)))
m = readvec('m.dat',shape=(len(y),len(x)))
Hexact = readvec('Hexact.dat',shape=(len(y),len(x)))
H = readvec('H.dat',shape=(len(y),len(x)))
residual = readvec('residual.dat',shape=(len(y),len(x)))

figdebug = False
def figsave(name):
    if figdebug:
        print '  showing %s ... close window to proceed ...' % name
        plt.show()  # debug
    else:
        plt.savefig(name,dpi=200,bbox_inches='tight')
        print '  figure file %s generated ...' % name

x = x/1000.0
y = y/1000.0

if len(x) == len(y):
    fsize = (12,9)
else:
    fsize = (int(19.0 * float(len(x)) / float(len(y))),14)

def gets(H,b):
    Hdraft = -(910.0/1028.0) * H
    icebase = np.maximum(b,Hdraft)
    s = H + icebase
    return s

if args.observed:
    Hexlabel = 'observed'
else:
    Hexlabel = 'exact'

half = (x > 0)
def extracthalf(f):
    if args.half:
        fn = f[n][half]
    else:
        fn = f[n][:]
    return fn

if args.profile:
    plt.figure(figsize=(7,5))
    n = len(y) / 2    # FIXME: make adjustable?
    xoriginal = x.copy()
    if args.half:
        x = x[half]
    bn = extracthalf(b)
    Hexactn = extracthalf(Hexact)
    plt.hold(True)
    # first plot bed and exact H
    if bn.max() - bn.min() > 0.0:  # only plot bed if it is not constant
        if args.exactbed:
            if not args.half:
                print "ERROR: option --exactbed only makes sense with --half"
                sys.exit(9)
            xfine = np.linspace(0.0,x.max()*1000.0,2001)
            bbed, Hbed = bedstep_exact(xfine)
            plt.plot(xfine/1000.0,bbed,'k',lw=2.5)
            plt.plot(xfine/1000.0,gets(Hbed,bbed),'-k',label=Hexlabel,lw=0.8)
        else:
            plt.plot(x,bn,'k',lw=1.5)
            plt.plot(x,gets(Hexactn,bn),'-k',label=Hexlabel,lw=0.8)
    else:
        if args.exactdome:
            if not args.half:
                print "ERROR: option --exactdome only makes sense with --half"
                sys.exit(8)
            xfine = np.linspace(0.0,x.max()*1000.0,2001)
            Hdome = dome_exact(xfine)
            plt.plot(xfine/1000.0,Hdome,'-k',label=Hexlabel,lw=0.8)
            if args.blowup:
                xfineblowup = xfine[(xfine >= 700.0e3) & (xfine <= 800.0e3)]
                Hdomeblowup = Hdome[(xfine >= 700.0e3) & (xfine <= 800.0e3)]
        else:
            plt.plot(x,Hexactn,'-k',label=Hexlabel,lw=0.8)
    # now plot extra H if any
    extrastylelist = ['+k', 'xk']
    if extraH:
        mainlabel = Hlabels[0]
    else:
        mainlabel = 'M*'
    if extraH and not args.blowup:  # don't plot extra_H in main figure ... see blowup below
        for j in range(len(extrastylelist)):  # FIXME ignore extra_H beyond 2
            nextH = readvec(Hlist[j],shape=(len(y),len(xoriginal)))
            nextHn = extracthalf(nextH)
            plt.plot(x,gets(nextHn,bn),extrastylelist[j],label=Hlabels[j+1],markersize=[10.0,8.0][j])
    # plot main numerical H last
    Hn = extracthalf(H)
    if bn.max() - bn.min() > 0.0:  # only plot bed if it is not constant
        plt.plot(x,gets(Hn,bn),'ok',label=mainlabel,markersize=6.0)
    else:
        plt.plot(x,Hn,'ok',label=mainlabel,markersize=6.0)
    # finish up with labels etc.
    plt.hold(False)
    plt.xlabel('x  (km)')
    plt.ylabel('z  (m)')
    plt.grid(True)
    plt.legend(fontsize=12.0)
    # fill in blowup
    if args.blowup:
        if not args.exactdome:
            print "ERROR: option --blowup only valid if also --exactdome"
            sys.exit(7)
        # put box around margin in original axes, and draw lines
        plt.hold(True)
        plt.plot([700.0, 700.0, 800.0, 800.0],[0.0, 1000.0, 1000.0, 0.0],'-k',lw=1.0)
        plt.plot([275.0, 685.0],[2390.0, 1025.0],'-k',lw=0.4)
        plt.plot([275.0, 620.0],[150.0, 25.0],'-k',lw=0.4)
        plt.hold(False)
        # add axes and replot quantities
        a = plt.axes([0.15, 0.13, 0.2, 0.45], axisbg='w')
        plt.plot(xfineblowup/1000.0,Hdomeblowup,'-k',label=Hexlabel,lw=0.8)
        plt.hold(True)
        iblowup = (x >= 700.0) & (x <= 800.0)
        xblowup = x[iblowup]
        if extraH:
            nextx = readvec('xfiner.dat')
            nexty = readvec('yfiner.dat')
            nextH = readvec('Hfiner.dat',shape=(len(nexty),len(nextx)))
            nextx /= 1000.0
            inext = (nextx >= 700.0) & (nextx <= 800.0)
            nextn = len(nexty) / 2
            plt.plot(nextx[inext],nextH[nextn][inext],'dk',markersize=5.0)
        plt.plot(xblowup,Hn[iblowup],'ok',markersize=7.0)
        plt.hold(False)
        plt.axis([700.0, 800.0, -100.0, 1000.0])
        plt.setp(a, xticks=[], yticks=[])  # no ticks on inset
    # output it
    figsave('HHexact1d.pdf')

if args.map:
    plt.figure(figsize=fsize)
    plt.pcolormesh(x,y,H)
    plt.axis('tight')
    plt.colorbar()
    plt.title('thickness solution H (m) with min=%.2f, max=%.2f' % (H.min(),H.max()))
    figsave('H.png')

    plt.figure(figsize=fsize)
    plt.pcolormesh(x,y,b)
    plt.axis('tight')
    if (b.max() > b.min()):
        plt.colorbar()
    plt.title('bed elevation b (m) with min=%.2f, max=%.2f' % (b.min(),b.max()))
    figsave('b.png')

    plt.figure(figsize=fsize)
    s = gets(H,b)
    plt.pcolormesh(x,y,s)
    plt.axis('tight')
    plt.colorbar()
    plt.title('surface elevation s (m) with min=%.2f, max=%.2f' % (s.min(),s.max()))
    figsave('s.png')

    plt.figure(figsize=fsize)
    m = m * 31556926.0
    plt.pcolormesh(x,y,m)
    plt.axis('tight')
    plt.colorbar()
    plt.title('surface mass balance m (m/a) with min=%.2f, max=%.2f' % (m.min(),m.max()))
    figsave('m.png')

    plt.figure(figsize=fsize)
    Herror = H - Hexact
    plt.pcolormesh(x,y,Herror)
    plt.axis('tight')
    plt.colorbar()
    plt.title('thickness error H-Hexact (m) with min=%.2f, max=%.2f' % (Herror.min(),Herror.max()))
    figsave('Herror.png')

    plt.figure(figsize=fsize)
    plt.pcolormesh(x,y,np.minimum(residual,0.0))
    plt.axis('tight')
    plt.colorbar()
    plt.title('VI-chopped residual r with min=%.2e so only %.2e <= r <= 0 shown' % \
              (residual.min(),residual.min()))
    figsave('residual.png')

