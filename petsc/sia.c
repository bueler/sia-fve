/* (C) 2015 Ed Bueler */

#include "sia.h"

// ******** FUNCTIONS ********

/*
Compute delta, a factor of both the diffusivity D and the pseudo-velocity W,
from values of thickness and surface gradient:
   delta = Gamma |grad H + grad b|^{n-1}
Also applies power-regularization part of continuation scheme.
*/
PetscReal getdelta(Grad gH, Grad gb, PetscReal Gamma, PetscReal n) {
    const PetscReal sx  = gH.x + gb.x,
                    sy  = gH.y + gb.y;
    return Gamma * PetscPowReal(sx*sx + sy*sy,(n-1.0)/2);
}

Grad getW(PetscReal delta, Grad gb) {
    Grad W;
    W.x = - delta * gb.x;
    W.y = - delta * gb.y;
    return W;
}

/* In continuation scheme, n(1)=1.0 and n(0)=n. */
PetscReal ncont(const AppCtx *user) {
  return (1.0 - user->eps) * user->n + user->eps * 1.0;
}

/* In continuation scheme, D(1)=D_0 and D(0)=D. */
PetscReal Dcont(PetscReal delta, PetscReal H, const AppCtx *user) {
  const PetscReal n = ncont(user);
  return (1.0-user->eps) * delta * PetscPowReal(H,n+2.0) + user->eps * user->D0;
}

/* See sia.h for doc. */
PetscReal getflux(Grad gH, Grad gb, PetscReal H, PetscReal Hup,
                  PetscBool xdir, AppCtx *user) {
  const PetscReal n     = ncont(user),
                  delta = getdelta(gH,gb,user->Gamma,n),
                  D     = Dcont(delta,H,user);
  const Grad      W     = getW(delta,gb);
  user->maxD = PetscMax(user->maxD, D);
  if (xdir)
      return - D * gH.x + W.x * PetscPowReal(Hup,n+2.0);
  else
      return - D * gH.y + W.y * PetscPowReal(Hup,n+2.0);
}


// ******** DERIVATIVES ********

/*
Regularized derivative of delta with respect to nodal value H_l, at point (x,y)
on element u,v:
   d delta / dl = (d delta / d gH.x) * (d gH.x / dl) + (d delta / d gH.y) * (d gH.y / dl),
but
   d delta / d gH.x = Gamma * (n-1) * |gH + gb|^{n-3.0} * (gH.x + gb.x)
   d delta / d gH.y = Gamma * (n-1) * |gH + gb|^{n-3.0} * (gH.y + gb.y)
However, power n-3.0 can be negative, so generating NaN in areas where the
surface gradient gH+gb is zero or tiny.
*/
PetscReal DdeltaDl(Grad gH, Grad gb, Grad dgHdl, const AppCtx *user) {
    const PetscReal sx  = gH.x + gb.x,
                    sy  = gH.y + gb.y,
                    n   = ncont(user),
                    slopesqr = sx * sx + sy * sy + user->slopeeps * user->slopeeps,
                    tmp = user->Gamma * (n-1) * PetscPowReal(slopesqr,(n-3.0)/2);
    return tmp * sx * dgHdl.x + tmp * sy * dgHdl.y;
}

/*
Derivative of continuation diffusivity D with respect to nodal value H_l, at
point (x,y) on element u,v.  Since  Dcont = (1-eps) delta H^{n+2} + eps D0:
   d Dcont / dl = (1-eps) [ (d delta / dl) H^{n+2} + delta (n+2) H^{n+1} (d H / dl) ]
                = (1-eps) H^{n+1} [ (d delta / dl) H + delta (n+2) (d H / dl) ]
*/
PetscReal DDcontDl(PetscReal delta, PetscReal ddeltadl, PetscReal H, PetscReal dHdl,
                   const AppCtx *user) {
    const PetscReal n = ncont(user);
    const PetscReal Hpow = PetscPowReal(H,n+1.0);
    return (1.0-user->eps) * Hpow * ( ddeltadl * H + delta * (n+2.0) * dHdl );
}

/*
Derivative of pseudo-velocity W with respect to nodal value H_l, at point
(x,y) on element u,v.  Since  W = - delta * grad b:
   d W.x / dl = - (d delta / dl) * gb.x
   d W.y / dl = - (d delta / dl) * gb.y
*/
Grad DWDl(PetscReal ddeltadl, Grad gb) {
    Grad dWdl;
    dWdl.x = - ddeltadl * gb.x;
    dWdl.y = - ddeltadl * gb.y;
    return dWdl;
}

/* See sia.h for doc. */
PetscReal DfluxDl(Grad gH, Grad gb, Grad dgHdl,
                  PetscReal H, PetscReal dHdl, PetscReal Hup, PetscReal dHupdl,
                  PetscBool xdir, const AppCtx *user) {
    const PetscReal n        = ncont(user),
                    delta    = getdelta(gH,gb,user->Gamma,n),
                    ddeltadl = DdeltaDl(gH,gb,dgHdl,user);
    const Grad      W        = getW(delta,gb),
                    dWdl     = DWDl(ddeltadl,gb);
    const PetscReal D        = Dcont(delta,H,user),
                    dDdl     = DDcontDl(delta,ddeltadl,H,dHdl,user),
                    Huppow   = PetscPowReal(Hup,n+1.0),
                    Huppow2  = Huppow * Hup,
                    dHuppow  = (n+2.0) * Huppow * dHupdl;
    if (xdir)
        return - dDdl * gH.x - D * dgHdl.x + dWdl.x * Huppow2 + W.x * dHuppow;
    else
        return - dDdl * gH.y - D * dgHdl.y + dWdl.y * Huppow2 + W.y * dHuppow;
}

