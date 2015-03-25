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

/*
int ddeltadgH(Grad gH, Grad gb, PetscReal Gamma, PetscReal n, PetscReal *ddeltadH1, PetscReal *ddeltadH2) {
FIXME
    const PetscReal sx  = gH.x + gb.x,
                    sy  = gH.y + gb.y;
    return Gamma * PetscPowReal(sx*sx + sy*sy,(n-1.0)/2);
}
*/

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

