/* (C) 2015 Ed Bueler */

/* Implementation of fundamental operations on bilinear (Q1 FEM) interpolants. */

#ifndef Q1OP_H_
#define Q1OP_H_

#include <petscvec.h>
#include "mahaffyctx.h"  // for Grad and AppCtx

// evaluate scalar field f(x,y)  at point (x,y), on element j,k, where
// locx = x-x_j, locy = y-y_k, using DMDA Vec f[k][j]
PetscErrorCode fieldatpt(PetscInt j, PetscInt k, PetscReal locx, PetscReal locy,
                         PetscReal **f, const AppCtx *user, PetscReal *f_atpt);

// evaluate gradient  (\nabla f)(x,y)  at point (x,y), on element j,k, where
// locx = x-x_j, locy = y-y_k, using DMDA Vec f[k][j]
PetscErrorCode gradatpt(PetscInt j, PetscInt k, PetscReal locx, PetscReal locy,
                        PetscReal **f, const AppCtx *user, Grad *gradf_atpt);

#endif // Q1OP_H_
