/* (C) 2015 Ed Bueler */

/* Implementation of SIA flux calculation. */

#ifndef SIA_H_
#define SIA_H_

#include <petsc.h>
#include "mahaffyctx.h"  // for Grad and AppCtx

/*
Evaluate a component of the flux from gradients and thickness.  Note:
   D = delta * H^{n+2}      (positive scalar)
   W = - delta * grad b     (vector)
so:
   q = D grad H + W H^{n+2} (vector)

If xdir==PETSC_TRUE then we compute
   q.x = - D gH.x + W.x H^{n+2}
but with upwinding, so we actually compute
   q.x = - D gH.x + W.x Hup^{n+2}
and similarly for y component if xdir==PETSC_FALSE.

This method also applies diffusivity- and power-regularization parts of the
continuation scheme.

Side effect: we update maximum of diffusivity (user->maxD).
*/
PetscReal getflux(Grad gH, Grad gb, PetscReal H, PetscReal Hup,
                  PetscBool xdir, AppCtx *user);

#endif // SIA_H_
