/* (C) 2015 Ed Bueler */

/* Implementation of fundamental operations on bilinear (Q1 FEM) interpolants. */

#ifndef Q1OP_H_
#define Q1OP_H_

#include <petscvec.h>
#include "mahaffyctx.h"  // for Grad only

/*
In the following functions
   * f(x,y) is represented by PetscReal** array for DMDA Vec, i.e. by {f[k][j]}
   * we are evaluating on element \square_{u,v}
   * use local coordinates:  xi = (x-x_u)/dx, eta = (y-y_v)/dy
       so (xi,eta) is in [0,1]x[0,1]
   * l indexes corner in counter-clockwise order:
       l=0 is (x_u,    y_v)
       l=1 is (x_{u+1},y_v)
       l=2 is (x_{u+1},y_{v+1})
       l=3 is (x_u,    y_{v+1})
*/

// evaluate scalar field  f(x,y)  at point (x,y) on element u,v
PetscReal fieldatpt(PetscInt u, PetscInt v, PetscReal xi, PetscReal eta, PetscReal **f);

// evaluate gradient  (\nabla f)(x,y)  at point (x,y) on element u,v
Grad      gradfatpt(PetscInt u, PetscInt v, PetscReal xi, PetscReal eta, PetscReal dx, PetscReal dy, PetscReal **f);

// evaluate derivative of scalar field  f(x,y)  with respect to nodal value f_l
// at point (x,y) on element u,v
PetscReal dfieldatpt(PetscInt l, PetscInt u, PetscInt v, PetscReal xi, PetscReal eta);

// evaluate derivative of gradient   (\nabla f)(x,y)  with respect to nodal value f_l
// at point (x,y) on element u,v
Grad      dgradfatpt(PetscInt l, PetscInt u, PetscInt v, PetscReal xi, PetscReal eta, PetscReal dx, PetscReal dy);

#endif // Q1OP_H_
