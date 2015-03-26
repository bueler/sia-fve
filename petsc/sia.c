/* (C) 2015 Ed Bueler */

#include "q1op.h"
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
  if (xdir == PETSC_TRUE) {
      return - D * gH.x + W.x * PetscPowReal(Hup,n+2.0);
  } else {
      return - D * gH.y + W.y * PetscPowReal(Hup,n+2.0);
  }
}


// ******** DERIVATIVES ********

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
Derivative of delta with respect to nodal value H_l, at point (x,y) on
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
Derivative of continuation diffusivity D with respect to nodal value H_l, at
point (x,y) on element u,v.  Since  Dcont = (1-eps) delta H^{n+2} + eps D0:
   d Dcont / dl = (1-eps) [ (d delta / dl) H^{n+2} + delta (n+2) H^{n+1} (d H / dl) ]
                = (1-eps) H^{n+1} [ (d delta / dl) H + delta (n+2) (d H / dl) ]
*/
PetscReal DDcontDl(PetscInt l,
                   PetscInt u, PetscInt v, PetscReal locx, PetscReal locy,
                   PetscReal delta,
                   PetscReal **H, PetscReal **b, const AppCtx *user) {
    const PetscReal n = ncont(user);
    // get d/dl values:
    const PetscReal Hatpt    = fieldatpt(u,v,locx,locy,H,user),
                    Hpow     = PetscPowReal(Hatpt,n+1.0),
                    dHdl     = dfieldatpt(l,u,v,locx,locy,user),
                    ddeltadl =   DdeltaDl(l,u,v,locx,locy,H,b,user);
    return (1.0-user->eps) * Hpow * ( ddeltadl * Hatpt + delta * (n+2.0) * dHdl );
}

/*
Derivative of pseudo-velocity W with respect to nodal value H_l, at point
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
// where D = Dcont and n = ncont
/*
PetscReal DfluxDl(PetscInt l,
                  PetscInt u, PetscInt v, PetscReal locx, PetscReal locy,
                  PetscReal Hup, PetscReal dHupdl,
                  PetscReal **H, PetscReal **b, const AppCtx *user) {
    const Grad gH = gradfatpt(u,v,locx,locy,H,user),
               gb = gradfatpt(u,v,locx,locy,b,user),
               dgHdl = dgradfatpt(l,u,v,locx,locy,user),
               dWdl  = DWDl(l,u,v,locx,locy,H,b,user);
    const PetscReal n     = ncont(user),
                    Hatpt = fieldatpt(u,v,locx,locy,H,user),
                    delta = getdelta(gH,gb,user->Gamma,n),
                    D     = Dcont(delta,Hatpt,user);
    const Grad W = getW(delta,gb,user);
    // get d/dl values:
    const PetscReal dDdl = DDcontDl(l,u,v,locx,locy,delta,H,b,user),
    dHdl     = dfieldatpt(l,u,v,locx,locy,user),
                    ddeltadl =   DdeltaDl(l,u,v,locx,locy,H,b,user);
    return (1.0-user->eps) * Hpow * ( ddeltadl * Hatpt + delta * (n+2.0) * dHdl );
}
*/
