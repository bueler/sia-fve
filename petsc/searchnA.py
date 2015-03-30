#!/usr/bin/env python
#
# (C) 2015 Ed Bueler
#
# Example of usage: FIXME

import argparse

parser = argparse.ArgumentParser(description='Optimize two-parameter functions by grid search or Nelder-Mead.')
parser.add_argument("--fprint", action="store_true", help="function will print at each evaluation")
parser.add_argument("--nm", action="store_true", help="use Nelder-Mead method")
parser.add_argument("--rosen", action="store_true", help="use rosen example")
parser.add_argument("--show", action="store_true", help="in grid search, plot surface")
args = parser.parse_args()

import numpy as np
import matplotlib.pyplot as plt

_status_message = {'success': 'Optimization terminated successfully.',
                   'maxiter': 'Maximum number of iterations has been '
                              'exceeded.'}

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
        #plt.pcolor(yy,xx[::-1],ff)
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
def neldermead(func, sim0, xtol=1e-4, ftol=1e-4, maxiter=None, disp=False):
    """
    Minimization of scalar function of one or more variables using the
    Nelder-Mead algorithm.

    Options for the Nelder-Mead algorithm are:
        xtol : float
            Relative error in solution `xopt` acceptable for convergence.
        ftol : float
            Relative error in ``fun(xopt)`` acceptable for convergence.
        maxiter : int
            Maximum number of iterations to perform.
        disp : bool
            Set to True to print convergence messages.
    """

    Np1 = sim0.shape[0]
    N = sim0.shape[1]
    if not N+1==Np1:
        raise ValueError("Initial simplex must have N+1 rows and N columns.")

    if maxiter is None:
        maxiter = N * 200

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
        simsize = np.max(np.ravel(np.abs(sim[1:] - sim[0])))
        fnorm = np.max(np.abs(fsim))
        #print simsize, fnorm
        if (simsize <= xtol and fnorm/fnorminit <= ftol):
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
        if disp:
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
        res, _, _, _ = neldermead(rosen, sim0, xtol=1e-1, ftol=1e-3, disp=True)
    else:
        res, _, _, _ = gridsearch(rosen, sim0, nx=3, ny=3)

else:
    import subprocess

    secpera = 31556926.0
    A0      = 1.0e-16/secpera

    def mahaffy(x):
        #cmd = "./mahaffy -mah_averr -mah_n %.4f -mah_A %.4e -mah_D0 10.0" % (x[0],x[1])
        cmd = "../mahaffy -mah_read grn10km.dat -mah_D0 40 -mah_Neps 10 -mah_averr -snes_type vinewtonssls -mah_n %.4f -mah_A %.4e" % (x[0],x[1])
        mahout = subprocess.check_output(cmd.split())
        if "DIV" in mahout:
            f = np.inf
            cmd = mahout + cmd
        else:
            f = float(mahout)
        if args.fprint:
            print "  f([%.4f %.4e]) = %.10f  from %s" % (x[0],x[1],f,cmd)
        return f

    sim0 = np.array([[3.0,A0],
                     [3.05,A0],
                     [3.0,1.15*A0]])

    if args.nm:
        res, _, _, _ = neldermead(mahaffy, sim0, xtol=1e-1, ftol=1e-3, disp=True)
    else:
        res, _, _, _ = gridsearch(mahaffy, sim0, nx=5, ny=5)

print res

