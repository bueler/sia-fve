# (C) 2015 Ed Bueler
#
# Compute exact solutions for use in profile figures.

import numpy as np

# exact thickness from Bueler (2003)
# compare DomeExactThickness() in exactsia.c
def dome_exact(x):
    L  = 750.0e3
    H0 = 3600.0
    n  = 3.0
    mm = 1.0 + 1.0 / n
    qq = n / (2.0 * n + 2.0)
    CC = H0 / (1.0 - 1.0 / n)**qq
    Hexact = 0.0 * x.copy()
    s   = x[x < L] / L
    tmp = mm * s - (1.0 / n) + (1.0 - s)**mm - s**mm
    Hexact[x < L] = CC * tmp**qq
    return Hexact

# bedrock elevation and exact thickness from Jarosch et al 2013
# formulas (56)--(59) in Jarosch et al (2013)
# compare BedStepExactThickness() in exactsia.c
def bedstep_exact(x):
    secpera = 3.1556926e7
    xm   = 20000.0
    xs   = 7000.0
    b0   = 500.0
    m0   = 2.0 / secpera
    n    = 3.0
    g    = 9.81
    rho  = 910.0
    A    = 1.0e-16/secpera
    rg   = rho * g
    ninv = 1.0 / n
    p    = n / (2.0*n+2.0)
    pinv = 1.0 / p
    q    = (2.0*n-1.0) / n
    den  = 6.0 * n * rg * (2.0*A)**ninv * xm**q
    Ch   = (2.0*n+2.0) * ((n+2.0)*m0)**ninv / den
    hsplus  = (Ch * (xm + 2.0*xs) * (xm - xs) * (xm - xs))**p
    hsminus = np.maximum(hsplus - b0, 0.0)
    Ca      = hsminus**pinv - hsplus**pinv
    # compute exact bed
    bexact = 0.0 * x.copy()
    bexact[abs(x)<xs] = 500.0
    # finally compute exact thickness
    Hexact = 0.0 * x.copy()
    absx = abs(x)
    tmp = Ch * (xm + 2.0 * absx) * (xm - absx) * (xm - absx)
    ileft  = (absx < xs)  & (absx < xm)
    iright = (absx >= xs) & (absx < xm)
    Hexact[ileft]  = (Ca + tmp[ileft])**p
    Hexact[iright] = tmp[iright]**p
    return bexact, Hexact

