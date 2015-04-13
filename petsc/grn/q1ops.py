# (C) 2015 Ed Bueler
# Utilities to do linear and binlinear interpolation.

import numpy as np

# coarse-to-fine linear interpolation in 1D:
#     c = coarse, f = fine, r = refinement factor
def lininterp(c,r):
    f = np.zeros(r * (len(c)-1) + 1)
    for j in range(len(c)-1):
        for s in range(r):
            lam = float(s)/float(r)
            f[r*j + s] = (1.0 - lam) * c[j] + lam * c[j+1]
    f[-1] = c[-1]
    return f

# coarse-to-fine bilinear interpolation plus padding in 2D
def quadinterp(c,r):
    f = np.zeros((r*(np.shape(c)[0]-1)+1,r*(np.shape(c)[1]-1)+1))
    for j in range(np.shape(c)[0]-1):
        for k in range(np.shape(c)[1]-1):
            for s in range(r):
                ls = float(s)/float(r)
                for t in range(r):
                    lt = float(t)/float(r)
                    f[r*j+s,r*k+t] =   (1.0-ls) * (1.0-lt) * c[j,k] \
                                     + ls       * (1.0-lt) * c[j+1,k] \
                                     + (1.0-ls) * lt       * c[j,k+1] \
                                     + ls       * lt       * c[j+1,k+1]
    f[:,-1] = lininterp(c[:,-1],r)
    f[-1,:] = lininterp(c[-1,:],r)
    return f

