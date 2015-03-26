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
PetscReal getflux(Grad gH, Grad gb,
                  PetscReal H,
                  PetscReal Hup,
                  PetscBool xdir, AppCtx *user);

/*
Derivative of flux (vector) with respect to nodal value H_l, at point
(x,y) on element u,v.  Since (above)
   q.x = - D gH.x + W.x Hup^{n+2},
we have:
   d q.x / dl = - (d D / dl) gH.x - D (d gH.x / dl)
                + (d W.x / dl) Hup^{n+2} + W.x (n+2) Hup^{n+1} (d Hup / dl),
where D = Dcont and n = ncont, and similar for d q.y / dl.

Calls hidden methods to compute  D,  W,  d D / dl,  d W / dl.
*/
PetscReal DfluxDl(Grad gH, Grad gb, Grad dgHdl,
                  PetscReal H, PetscReal dHdl,
                  PetscReal Hup, PetscReal dHupdl,
                  PetscBool xdir, const AppCtx *user);

#endif // SIA_H_
