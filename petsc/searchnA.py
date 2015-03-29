#!/usr/bin/env python
#
# (C) 2015 Ed Bueler
#
# Example of usage: FIXME

import argparse

parser = argparse.ArgumentParser(description='Optimize two-parameter functions by Nelder-Mead.')
parser.add_argument("--rosen", action="store_true", help="use rosen example")
args = parser.parse_args()

import sys
import numpy

_status_message = {'success': 'Optimization terminated successfully.',
                   'maxiter': 'Maximum number of iterations has been '
                              'exceeded.'}

# FIXME: add this
#def gridsearch(func, sim0, nx=3, ny=3, disp=False):

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
    fsim = numpy.zeros((N + 1,))
    for k in range(0, N+1):
        y = numpy.array(sim[k,:], copy=True)
        fsim[k] = func(y)

    fnorminit = numpy.max(numpy.abs(fsim))

    #print sim
    #print fsim

    rho = 1.0
    chi = 2.0
    psi = 0.5
    sigma = 0.5
    one2np1 = list(range(1, N + 1))

    # sort so sim[0,:] has the lowest function value
    ind = numpy.argsort(fsim)
    fsim = numpy.take(fsim, ind, 0)
    sim = numpy.take(sim, ind, 0)

    iterations = 1

    while (iterations < maxiter):
        simsize = numpy.max(numpy.ravel(numpy.abs(sim[1:] - sim[0])))
        fnorm = numpy.max(numpy.abs(fsim))
        #print simsize, fnorm
        if (simsize <= xtol and fnorm/fnorminit <= ftol):
            break

        xbar = numpy.add.reduce(sim[:-1], 0) / N
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

        ind = numpy.argsort(fsim)
        sim = numpy.take(sim, ind, 0)
        fsim = numpy.take(fsim, ind, 0)
        iterations += 1

    x = sim[0]
    fval = numpy.min(fsim)
    warnflag = 0

    if iterations >= maxiter:
        warnflag = 2
        msg = _status_message['maxiter']
        if disp:
            print('Warning: ' + msg)
    else:
        if disp:
            print(_status_message['success'])
            print("         Current function value: %f" % fval)
            print("         Iterations: %d" % iterations)
            print sim

    return x


if args.rosen:
    def rosen(x):
        """The Rosenbrock function"""
        f = sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)
        print "  f([%12.9f %12.9f]) = %8e" % (x[0],x[1],f)
        return f
    sim0 = numpy.array([[1.3, 0.7],
                        [1.4, 0.7],
                        [1.3, 0.8]])
    res = neldermead(rosen, sim0, xtol=1e-1, ftol=1e-3, disp=True)

else:
    secpera = 31556926.0
    A0      = 1.0e-16/secpera

    raise NotImplementedError("Mahaffy case not implemented.")

    # FIXME: build function by calling mahaffy with options
    #   -mah_n n -mah_A A
    # and returning L^1 error norm as function value
    def mahaffy(x):
        #sys.system("./mahaffy -mah_n ?? -mah_A ?? ...")
        # FIXME extract L^1 norm as f
        print "  f([%12.9f, %12.9f A0]) = %8e" % (x[0],x[1]/A0,f)
        return f
    sim0 = numpy.array([[3.0,A0],
                        [3.5,A0],
                        [3.0,2.0*A0]])
    res = neldermead(mahaffy, sim0, xtol=1e-1, ftol=1e-3, disp=True)

print res

