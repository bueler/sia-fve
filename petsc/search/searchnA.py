#!/usr/bin/env python
#
# (C) 2015 Ed Bueler

import argparse
import subprocess
import numpy as np
from optimizers import neldermead2, gridsearch
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=
'''
Optimize two-parameter objective functions by grid search (default) or
Nelder-Mead.  Ddesigned to run ../mahaffy with various values of
n and A and use the numerical error output as the objective value.
''')

parser.add_argument("--fprint", action="store_true",
                    help="function will print at each evaluation")

ngroup = parser.add_argument_group("Nelder-Mead method")
ngroup.add_argument("--disp", action="store_true",
                    help="show simplex and convergence condition")
ngroup.add_argument("--nm", action="store_true",
                    help="use Nelder-Mead method")

ggroup = parser.add_argument_group("grid search method")
ggroup.add_argument("--nx", type=int, default=3,
                    help="use  2 nx + 1  points in x direction (= n direction; default: %(default)d)")
ggroup.add_argument("--ny", type=int, default=3,
                    help="use  2 ny + 1  points in y direction (= A direction; default: %(default)d)")
ggroup.add_argument("--show", action="store_true",
                    help="color plot of objective evaluations")

rgroup = parser.add_argument_group("test objective")
rgroup.add_argument("--rosen", action="store_true",
                    help="use Rosenbrock example objective function for testing")

mgroup = parser.add_argument_group("mahaffy output objective")
mgroup.add_argument("--mah", metavar = 'STR',
                    default="-cs_D0 40 -cs_end 10 -snes_type vinewtonssls -pc_type asm -sub_pc_type lu",
                    help="arguments to mahaffy (default: '%(default)s')")
mgroup.add_argument("--maxerr", action="store_true",
                    help="use -mah_maxerr instead of -mah_averr")
mgroup.add_argument("--more", default="", metavar = 'STR',
                    help="additional arguments to mahaffy")
mgroup.add_argument("--mpi", type=int, default=1, metavar='N',
                    help="number of processes in 'mpiexec -n N' (default: %(default)d)")
mgroup.add_argument("--read", metavar='FILE', default='',
                    help="file to read using -mah_read")
mgroup.add_argument("--restart", action="store_true",
                    help="use -mah_dump foo/ and -mah_read foo/unnamed.dat -mah_readinitial")

args = parser.parse_args()

if args.nm:
    method = "Nelder-Mead"
    if args.show:
        raise ValueError("invalid option --show with Nelder-Mead method")
else:
    method = "grid search"
    if args.disp:
        raise ValueError("invalid option --disp with grid search")

print "using %s to optimize output from" % method,

if args.rosen:
    print "Rosenbrock test function ..."
    if len(args.more)>0 or args.read or (args.mpi>1):
        raise ValueError("invalid option when optimizing Rosenbrock test function")
else:
    # conversion of options into mahaffy command
    mahargs = args.mah
    if len(args.more) > 0:
        mahargs = mahargs + " " + args.more
    if args.maxerr:
        maherr = "-mah_maxerr"
    else:
        maherr = "-mah_averr"
    if args.mpi > 1:
        mahmpistr = "mpiexec -n %d " % args.mpi
    else:
        mahmpistr = ""
    cmdroot = "%s../mahaffy -mah_read %s %s %s " % (mahmpistr, args.read, maherr, mahargs)
    cmdfmt = cmdroot + "-mah_n N -mah_A A"
    print "mahaffy command:"
    print "    f([N A]) = { " + cmdfmt + " }"

if args.rosen:
    def rosen(x):
        """The Rosenbrock function"""
        f = sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)
        if args.fprint:
            print "  f([%12.9f %12.9f]) = %8e" % (x[0],x[1],f)
        return f
    myfunc = rosen
    sim0 = np.array([[1.05, 0.93],
                     [1.15, 0.93],
                     [1.05, 1.0]])
    myxtol=[5e-2,5e-2]
    myftol=1e-3
else:
    secpera = 31556926.0
    A0      = 1.0e-16/secpera
    def mahaffy(x):
        if args.fprint:
            print "f([%.4f %.4e]) = " % (x[0],x[1]),
        cmd = cmdroot + "-mah_n %.4f -mah_A %.4e" % (x[0], x[1])
        if args.restart:
            import os.path
            if not os.path.exists('foo/'):
                subprocess.call("mkdir -p foo/".split())
            cmd += " -mah_dump foo/"
            if os.path.exists('foo/unnamed.dat'):
                cmd += " -mah_read foo/unnamed.dat -mah_readinitial -cs_start 8"
        mahout = subprocess.check_output(cmd.split())
        if "DIV" in mahout:
            f = np.inf
            if args.fprint:
                print "inf  ... because returned %s" % mahout
        else:
            f = float(mahout)
            if args.fprint:
                print "%.10f" % f
        return f
    myfunc = mahaffy
    sim0 = np.array([[3.0,A0],
                     [3.05,A0],
                     [3.0,1.15*A0]])
    myxtol=[0.01, 0.1e-24]
    myftol=np.inf

if args.nm:
    res, _, _, _ = neldermead2(myfunc, sim0, xtol=myxtol, ftol=myftol, disp=args.disp)
else:
    nx = args.nx
    ny = args.ny
    dx = sim0[1,0] - sim0[0,0]
    dy = sim0[2,1] - sim0[0,1]
    print "running on grid of %d by %d points:" % (2*nx+1,2*ny+1)
    print "    [%.2e, ..., %.2e] x [%.2e, ..., %.2e]" % \
       (sim0[0,0]-dx*nx, sim0[0,0]+dx*nx, sim0[0,1]-dy*ny, sim0[0,1]+dy*ny)
    res, _, _, _ = gridsearch(myfunc, sim0, nx=args.nx, ny=args.ny, show=args.show)
    if args.show:
        if not args.rosen:
            plt.xlabel('ice softness A')
            plt.ylabel('Glen exponent n')
        plt.show()

print res

