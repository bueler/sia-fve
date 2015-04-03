# (C) 2015 Ed Bueler

import numpy as np
import matplotlib.pyplot as plt

def gridsearch(func, sim0, nx=3, ny=3, show=False):
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

    if show:
        plt.imshow(ff, interpolation='nearest', origin='lower',
                   extent=[yy.min()-dy/2,yy.max()+dy/2,xx.min()-dx/2,xx.max()+dx/2])
        plt.colorbar()
        plt.axis('tight')

    # report results like neldermead()
    return xymin, iterations, funcmin, 0

# following is modified _minimize_neldermead() from
#   https://github.com/scipy/scipy/blob/v0.14.0/scipy/optimize/optimize.py
def neldermead2(func, sim0, xtol=[1e-4,1e-4], ftol=1e-4, maxiter=np.inf, disp=False):
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
        if disp:
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
        print('Warning: Maximum number of iterations has been exceeded.')

    return x, iterations, fval, warnflag

