#!/usr/bin/env python
#
# (C) 2015 Ed Bueler

import argparse

parser = argparse.ArgumentParser(description='Optimize two-parameter functions of n and A by grid search (default) or Nelder-Mead.')
parser.add_argument("file", metavar='FILE', help="input file to -mah_read")
parser.add_argument("--disp", action="store_true", help="optimization method prints messages")
parser.add_argument("--fprint", action="store_true", help="function will print at each evaluation")
parser.add_argument("--mah", default="-mah_D0 40 -mah_Neps 10 -snes_type vinewtonssls",
                    help="arguments to mahaffy (default: '%(default)s')")
parser.add_argument("--maxerr", action="store_true", help="use -mah_maxerr instead of -mah_averr")
parser.add_argument("--more", default="", help="arguments to mahaffy to add")
parser.add_argument("--nm", action="store_true", help="use Nelder-Mead method")
parser.add_argument("--nx", type=int, default=3,
                    help="use  2 nx + 1  points in N direction in grid search (default: %(default)d)")
parser.add_argument("--ny", type=int, default=3,
                    help="use  2 ny + 1  points in A direction in grid search (default: %(default)d)")
parser.add_argument("--rosen", action="store_true", help="use Rosenbrock example objective function (for testing)")
parser.add_argument("--show", action="store_true", help="in grid search, plot surface")
args = parser.parse_args()

import numpy as np
import matplotlib.pyplot as plt

if args.nm:
    method = "Nelder-Mead"
else:
    method = "grid search"
mahargs = args.mah + args.more
if args.maxerr:
    maherr = "-mah_maxerr"
else:
    maherr = "-mah_averr"
cmdroot = "./mahaffy -mah_read %s %s %s " % (args.file, maherr, mahargs)
cmdfmt = cmdroot + "-mah_n N -mah_A A"

print "using %s to optimize output from" % method,
if args.rosen:
    print "Rosenbrock function ..."
else:
    print "command:"
    print "    f([N A]) = " + cmdfmt

_status_message = {'success': 'Optimization terminated successfully.',
                   'maxiter': 'Maximum number of iterations has been exceeded.'}

def gridsearch(func, sim0, nx=3, ny=3):
    Np1 = sim0.shape[0]
    N = sim0.shape[1]
    if not N+1==Np1:
        raise ValueError("Initial simplex must have N+1 rows and N columns.")
    if not N==2:
        raise ValueError("Only works for N==2.")

    # assume "simplex" set up a certain way: [[x,y], [x+dx,y], [x,y+dy]]
    dx = sim0[1,0] - sim0[0,0]
    dy = sim0[2,1] - sim0[0,1]

    # run grid search
    iterations = (2*nx+1) * (2*ny+1)
    ff = np.zeros((2*nx+1,2*ny+1))
    xx = np.linspace(sim0[0,0] - dx * nx, sim0[0,0] + dx * nx, 2*nx+1)
    yy = np.linspace(sim0[0,1] - dy * ny, sim0[0,1] + dy * ny, 2*ny+1)
    funcmin = np.inf
    xymin = [np.nan,np.nan]
    for j in range(2*nx+1):
        for k in range(2*ny+1):
            xy = np.array([xx[j],yy[k]])
            ff[j,k] = func(xy)
            if ff[j,k] <= funcmin:
                funcmin, xymin = ff[j,k], xy

    if args.show:
        plt.imshow(ff, interpolation='nearest', origin='lower',
                   extent=[yy.min()-dy/2,yy.max()+dy/2,xx.min()-dx/2,xx.max()+dx/2])
        plt.xlabel('ice softness A')
        plt.ylabel('Glen exponent n')
        plt.colorbar()
        plt.axis('tight')
        plt.show()

    # report results like neldermead()
    return xymin, iterations, funcmin, 0

# following is modified _minimize_neldermead() from
#   https://github.com/scipy/scipy/blob/v0.14.0/scipy/optimize/optimize.py
def neldermead2(func, sim0, xtol=[1e-4,1e-4], ftol=1e-4, maxiter=np.inf):
    """
    Minimization of scalar function of TWO variables using the Nelder-Mead algorithm.

    Options for the Nelder-Mead algorithm are:
        xtol : array
            Simplex dimensions (differences of x values) acceptable for convergence.
        rtol : float
            Relative reduction in func(x) acceptable for convergence.
        maxiter : int
            Maximum number of iterations to perform.
    """

    Np1 = sim0.shape[0]
    N = sim0.shape[1]
    if not N+1==Np1:
        raise ValueError("Initial simplex must have N+1 rows and N columns.")

    # evaluate f at pts of initial simplex
    sim = sim0.copy()
    fsim = np.zeros((N + 1,))
    for k in range(0, N+1):
        y = np.array(sim[k,:], copy=True)
        fsim[k] = func(y)

    fnorminit = np.max(np.abs(fsim))

    rho = 1.0
    chi = 2.0
    psi = 0.5
    sigma = 0.5
    one2np1 = list(range(1, N + 1))

    # sort so sim[0,:] has the lowest function value
    ind = np.argsort(fsim)
    fsim = np.take(fsim, ind, 0)
    sim = np.take(sim, ind, 0)

    iterations = 1

    while (iterations < maxiter):
        simdiffmax = np.abs(sim[1:,:] - sim[0,:]).max(0)
        fnorm = np.max(np.abs(fsim))
        if args.disp:
            print sim
            print "  ?: ", simdiffmax, " <= ", xtol, " AND ", fnorm/fnorminit, " <= ", ftol
        if (np.all(simdiffmax <= xtol) and fnorm/fnorminit <= ftol):
            break

        xbar = np.add.reduce(sim[:-1], 0) / N
        xr = (1 + rho) * xbar - rho * sim[-1]
        fxr = func(xr)
        doshrink = 0

        if fxr < fsim[0]:
            xe = (1 + rho * chi) * xbar - rho * chi * sim[-1]
            fxe = func(xe)

            if fxe < fxr:
                sim[-1] = xe
                fsim[-1] = fxe
            else:
                sim[-1] = xr
                fsim[-1] = fxr
        else:  # fsim[0] <= fxr
            if fxr < fsim[-2]:
                sim[-1] = xr
                fsim[-1] = fxr
            else:  # fxr >= fsim[-2]
                # Perform contraction
                if fxr < fsim[-1]:
                    xc = (1 + psi * rho) * xbar - psi * rho * sim[-1]
                    fxc = func(xc)

                    if fxc <= fxr:
                        sim[-1] = xc
                        fsim[-1] = fxc
                    else:
                        doshrink = 1
                else:
                    # Perform an inside contraction
                    xcc = (1 - psi) * xbar + psi * sim[-1]
                    fxcc = func(xcc)

                    if fxcc < fsim[-1]:
                        sim[-1] = xcc
                        fsim[-1] = fxcc
                    else:
                        doshrink = 1

                if doshrink:
                    for j in one2np1:
                        sim[j] = sim[0] + sigma * (sim[j] - sim[0])
                        fsim[j] = func(sim[j])

        ind = np.argsort(fsim)
        sim = np.take(sim, ind, 0)
        fsim = np.take(fsim, ind, 0)
        iterations += 1

    x = sim[0]
    fval = np.min(fsim)
    warnflag = 0

    if iterations >= maxiter:
        warnflag = 2
        msg = _status_message['maxiter']
        print('Warning: ' + msg)

    return x, iterations, fval, warnflag


if args.rosen:

    def rosen(x):
        """The Rosenbrock function"""
        f = sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)
        if args.fprint:
            print "  f([%12.9f %12.9f]) = %8e" % (x[0],x[1],f)
        return f

    sim0 = np.array([[1.3, 0.7],
                     [1.4, 0.7],
                     [1.3, 0.8]])
    if args.nm:
        res, _, _, _ = neldermead2(rosen, sim0, xtol=[1e-1,1e-1], ftol=1e-2)
    else:
        res, _, _, _ = gridsearch(rosen, sim0, nx=3, ny=3)

else:
    import subprocess
    secpera = 31556926.0
    A0      = 1.0e-16/secpera

    def mahaffy(x):
        cmd = cmdroot + "-mah_n %.4f -mah_A %.4e" % (x[0], x[1])
        mahout = subprocess.check_output(cmd.split())
        if "DIV" in mahout:
            f = np.inf
            if args.fprint:
                print "f([%.4f %.4e]) = inf  ... because returned %s" % (x[0],x[1],mahout)
        else:
            f = float(mahout)
            if args.fprint:
                print "f([%.4f %.4e]) = %.10f" % (x[0],x[1],f)
        return f

    sim0 = np.array([[3.0,A0],
                     [3.05,A0],
                     [3.0,1.15*A0]])
    if args.nm:
        res, _, _, _ = neldermead2(mahaffy, sim0, xtol=[0.1, 0.1e-24], ftol=np.inf)
    else:
        res, _, _, _ = gridsearch(mahaffy, sim0, nx=args.nx, ny=args.ny)

print res

