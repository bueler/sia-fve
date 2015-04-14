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

def locquadinterp(f, xi, eta):
    return   (1.0-xi) * (1.0-eta) * f[0,0]  + xi * (1.0-eta) * f[1,0] \
           + (1.0-xi) * eta       * f[0,1]  + xi       * eta * f[1,1]

# coarse-to-fine bilinear interpolation in 2D
def quadinterp(c,r):
    f = np.zeros((r*(np.shape(c)[0]-1)+1,r*(np.shape(c)[1]-1)+1))
    for j in range(np.shape(c)[0]-1):
        for k in range(np.shape(c)[1]-1):
            for s in range(r):
                xi = float(s)/float(r)
                for t in range(r):
                    eta = float(t)/float(r)
                    f[r*j+s,r*k+t] = locquadinterp(c[j:j+2,k:k+2],xi,eta)
    f[:,-1] = lininterp(c[:,-1],r)
    f[-1,:] = lininterp(c[-1,:],r)
    return f

