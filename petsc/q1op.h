/* (C) 2015 Ed Bueler */

/* Implementation of fundamental operations on bilinear (Q1 FEM) interpolants. */

#ifndef Q1OP_H_
#define Q1OP_H_

#include <petscvec.h>
#include "mahaffyctx.h"  // for Grad and AppCtx

/* In the following functions
   * f(x,y) is represented by DMDA Vec {f[k][j]}
   * use local coordinates:  locx = x-x_u, locy = y-y_v
   * l indexes corner in counter-clockwise order:
       l=0 is (x_u,    y_v)
       l=1 is (x_{u+1},y_v)
       l=2 is (x_{u+1},y_{v+1})
       l=3 is (x_u,    y_{v+1})
*/

// evaluate scalar field  f(x,y)  at point (x,y), on element u,v
PetscErrorCode fieldatpt(PetscInt u, PetscInt v, PetscReal locx, PetscReal locy,
                         PetscReal **f, const AppCtx *user, PetscReal *f_atpt);

// evaluate gradient  (\nabla f)(x,y)  at point (x,y), on element u,v
PetscErrorCode gradfatpt(PetscInt u, PetscInt v, PetscReal locx, PetscReal locy,
                         PetscReal **f, const AppCtx *user, Grad *gradf_atpt);

// evaluate derivative of scalar field  f(x,y)  with respect to nodal value f_l
// at point (x,y), on element u,v
PetscErrorCode dfieldatpt(PetscInt l,
                          PetscInt u, PetscInt v, PetscReal locx, PetscReal locy,
                          const AppCtx *user, PetscReal *dfdl_atpt);

// evaluate derivative of gradient   (\nabla f)(x,y)  with respect to nodal value f_l
// at point (x,y), on element u,v
PetscErrorCode dgradfatpt(PetscInt l,
                          PetscInt u, PetscInt v, PetscReal locx, PetscReal locy,
                          const AppCtx *user, Grad *dgradfdl_atpt);

#endif // Q1OP_H_
