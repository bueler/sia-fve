/* (C) 2015 Ed Bueler */

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


Grad ddeltadgH(Grad gH, Grad gb, PetscReal Gamma, PetscReal n) {
    const PetscReal sx  = gH.x + gb.x,
                    sy  = gH.y + gb.y,
                    tmp = Gamma * (n-1) * PetscPowReal(sx*sx + sy*sy,(n-3.0)/2);
    Grad ddH;
    ddH.x = tmp * sx;
    ddH.y = tmp * sy;
    return ddH;
}

// from above:
//   d delta / dl = ddH.x * (d gH.x / dl) + ddH.y * (d gH.y / dl)

/* FIXME
PetscErrorCode ddeltadl(PetscInt u, PetscInt v, PetscReal locx, PetscReal locy,
                         PetscReal **f, const AppCtx *user, PetscReal *f_atpt);
*/

// and since  D = delta H^{n+2}:
//   d D / dl = (d delta / dl) H^{n+2} + delta (n+2) H^{n+1} (d H / dl)
// and since  W = - delta * grad b:
//   d W.x / dl = - (d delta / dl) * b.x
//   d W.y / dl = - (d delta / dl) * b.y

// since  q.x = - D gH.x + W.x Hup^{n+2}  we have:
//   d q.x / dl = - (d D / dl) Hx - D (d gH.x / dl) + (d W.x / dl) Hup^{n+2} + W.x (n+2) Hup^{n+1} (d Hup / dl)

PetscReal getflux(Grad gH, Grad gb, PetscReal H, PetscReal Hup,
                  PetscBool xdir, AppCtx *user) {
  const PetscReal eps   = user->eps,
                  n     = (1.0 - eps) * user->n + eps * 1.0,
                  delta = getdelta(gH,gb,user->Gamma,n),
                  D     = (1.0-eps) * delta * PetscPowReal(H,n+2.0) + eps * user->D0;
  user->maxD = PetscMax(user->maxD, D);
  if (xdir == PETSC_TRUE) {
      const PetscReal  Wx = - delta * gb.x;
      return - D * gH.x + Wx * PetscPowReal(Hup,n+2.0);
  } else {
      const PetscReal  Wy = - delta * gb.y;
      return - D * gH.y + Wy * PetscPowReal(Hup,n+2.0);
  }
}

