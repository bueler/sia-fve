/* (C) 2015 Ed Bueler */

#include "q1op.h"
#include "sia.h"

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

/* In continuation scheme, n(1)=1.0 and n(0)=n. */
PetscReal ncont(const AppCtx *user) {
  return (1.0 - user->eps) * user->n + user->eps * 1.0;
}

/* In continuation scheme, D(1)=D_0 and D(0)=D. */
PetscReal Dcont(PetscReal delta, PetscReal H, const AppCtx *user) {
  const PetscReal n = ncont(user);
  return (1.0-user->eps) * delta * PetscPowReal(H,n+2.0) + user->eps * user->D0;
}

/*
Derivative of delta with respect to components of gH.  Returns Grad ddH with two
components:
   ddH.x = d delta / d gH.x
   ddH.y = d delta / d gH.y
*/
Grad DdeltaDgH(Grad gH, Grad gb, PetscReal Gamma, PetscReal n) {
    const PetscReal sx  = gH.x + gb.x,
                    sy  = gH.y + gb.y,
                    tmp = Gamma * (n-1) * PetscPowReal(sx*sx + sy*sy,(n-3.0)/2);
    Grad ddH;
    ddH.x = tmp * sx;
    ddH.y = tmp * sy;
    return ddH;
}

/*
Derivative of delta with respect to corner thickness H_l, at point (x,y) on
element u,v:
   d delta / dl = (d delta / d gH.x) * (d gH.x / dl) + (d delta / d gH.y) * (d gH.y / dl)
*/
PetscReal DdeltaDl(PetscInt l,
                   PetscInt u, PetscInt v, PetscReal locx, PetscReal locy,
                   PetscReal **H, PetscReal **b, const AppCtx *user) {
    const PetscReal n = ncont(user);
    const Grad gH    =    gradfatpt(u,v,locx,locy,H,user),
               gb    =    gradfatpt(u,v,locx,locy,b,user),
               dgHdl = dgradfatpt(l,u,v,locx,locy,user),
               ddH   = DdeltaDgH(gH,gb,user->Gamma,n);
    return ddH.x * dgHdl.x + ddH.y * dgHdl.y;
}

/*
Derivative of diffusivity D with respect to corner thickness H_l, at point (x,y) on
element u,v.  Since  D = delta H^{n+2}:
   d D / dl = (d delta / dl) H^{n+2} + delta (n+2) H^{n+1} (d H / dl)
*/
PetscReal DDDl(PetscInt l,
               PetscInt u, PetscInt v, PetscReal locx, PetscReal locy,
               PetscReal **H, PetscReal **b, const AppCtx *user) {
    const PetscReal n = ncont(user);
    // delta could be an argument to this method, in which case this not needed:
    const Grad gH = gradfatpt(u,v,locx,locy,H,user),
               gb = gradfatpt(u,v,locx,locy,b,user);
    const PetscReal delta = getdelta(gH,gb,user->Gamma,n);
    // get d/dl values:
    const PetscReal Hatpt    =    fieldatpt(u,v,locx,locy,H,user),
                    dHdl     = dfieldatpt(l,u,v,locx,locy,user),
                    ddeltadl =   DdeltaDl(l,u,v,locx,locy,H,b,user);
    return ddeltadl * PetscPowReal(Hatpt,n+2.0) + delta * (n+2.0) * PetscPowReal(Hatpt,n+1.0) * dHdl;
}

/*
Derivative of pseudo-velocity W with respect to corner thickness H_l, at point
(x,y) on element u,v.  Since  W = - delta * grad b:
   d W.x / dl = - (d delta / dl) * gb.x
   d W.y / dl = - (d delta / dl) * gb.y
*/
Grad DWDl(PetscInt l,
          PetscInt u, PetscInt v, PetscReal locx, PetscReal locy,
          PetscReal **H, PetscReal **b, const AppCtx *user) {
    const Grad gb = gradfatpt(u,v,locx,locy,b,user);
    const PetscReal ddeltadl = DdeltaDl(l,u,v,locx,locy,H,b,user);
    Grad dWdl;
    dWdl.x = - ddeltadl * gb.x;
    dWdl.y = - ddeltadl * gb.y;
    return dWdl;
}


//FIXME since  q.x = - D gH.x + W.x Hup^{n+2}  we have:
//   d q.x / dl = - (d D / dl) gH.x - D (d gH.x / dl) + (d W.x / dl) Hup^{n+2} + W.x (n+2) Hup^{n+1} (d Hup / dl)

/* See sia.h for doc. */
PetscReal getflux(Grad gH, Grad gb, PetscReal H, PetscReal Hup,
                  PetscBool xdir, AppCtx *user) {
  const PetscReal n     = ncont(user),
                  delta = getdelta(gH,gb,user->Gamma,n),
                  D     = Dcont(delta,H,user);
  user->maxD = PetscMax(user->maxD, D);
  if (xdir == PETSC_TRUE) {
      const PetscReal  Wx = - delta * gb.x;
      return - D * gH.x + Wx * PetscPowReal(Hup,n+2.0);
  } else {
      const PetscReal  Wy = - delta * gb.y;
      return - D * gH.y + Wy * PetscPowReal(Hup,n+2.0);
  }
}

