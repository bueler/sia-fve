/* (C) 2015 Ed Bueler */

#include "q1op.h"


int weights(PetscReal dx, PetscReal dy, PetscReal locx, PetscReal locy,
            PetscReal (*x)[4], PetscReal (*y)[4], PetscReal (*gx)[4], PetscReal (*gy)[4]) {
  const PetscReal xi  = locx / dx,
                  eta = locy / dy;
  if (x) {
    (*x)[0]  = 1.0 - xi;   (*x)[1]  = xi;         (*x)[2]  = xi;       (*x)[3]  = 1.0 - xi;
  }
  if (gx) {
    (*gx)[0] = - 1.0 / dx; (*gx)[1] = 1.0 / dx;   (*gx)[2] = 1.0 / dx; (*gx)[3] = - 1.0 / dx;
  }
  if (y) {
    (*y)[0]  = 1.0 - eta;  (*y)[1]  = 1.0 - eta;  (*y)[2]  = eta;      (*y)[3]  = eta;
  }
  if (gy) {
    (*gy)[0] = - 1.0 / dy; (*gy)[1] = - 1.0 / dy; (*gy)[2] = 1.0 / dy; (*gy)[3] = 1.0 / dy;
  }
  return 0;
}

PetscReal fieldatpt(PetscInt u, PetscInt v,         // (j,k) is the element (by lower-left corner)
                    PetscReal locx, PetscReal locy, // = (x,y) coords in element
                    PetscReal **f,                  // f[k][j] are node values
                    const AppCtx *user) {
  PetscReal x[4], y[4];
  weights(user->dx,user->dx,locx,locy,&x,&y,NULL,NULL);
  return x[0] * y[0] * f[v][u] + x[1] * y[1] * f[v][u+1] + x[2] * y[2] * f[v+1][u+1] + x[3] * y[3] * f[v+1][u];
}

Grad gradfatpt(PetscInt u, PetscInt v, PetscReal locx, PetscReal locy,
               PetscReal **f, const AppCtx *user) {
  PetscReal x[4], y[4], gx[4], gy[4];
  Grad gradf;
  weights(user->dx,user->dx,locx,locy,&x,&y,&gx,&gy);
  gradf.x =   gx[0] * y[0] * f[v][u]     + gx[1] * y[1] * f[v][u+1]
            + gx[2] * y[2] * f[v+1][u+1] + gx[3] * y[3] * f[v+1][u];
  gradf.y =    x[0] *gy[0] * f[v][u]     +  x[1] *gy[1] * f[v][u+1]
            +  x[2] *gy[2] * f[v+1][u+1] +  x[3] *gy[3] * f[v+1][u];
  return gradf;
}

PetscReal dfieldatpt(PetscInt l,
                     PetscInt u, PetscInt v, PetscReal locx, PetscReal locy,
                     const AppCtx *user) {
  PetscReal x[4], y[4];
  weights(user->dx,user->dx,locx,locy,&x,&y,NULL,NULL);
  return x[l] * y[l];
}

Grad dgradfatpt(PetscInt l,
                PetscInt u, PetscInt v, PetscReal locx, PetscReal locy,
                const AppCtx *user) {
  PetscReal x[4], y[4], gx[4], gy[4];
  Grad dgradfdl;
  weights(user->dx,user->dx,locx,locy,&x,&y,&gx,&gy);
  dgradfdl.x = gx[l] *  y[l];
  dgradfdl.y =  x[l] * gy[l];
  return dgradfdl;
}
