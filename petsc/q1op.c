/* (C) 2015 Ed Bueler */

#include "q1op.h"

#define WEIGHTS \
const PetscReal x[4] = { 1.0-xi,      xi,  xi, 1.0-xi}, \
                y[4] = {1.0-eta, 1.0-eta, eta,    eta};

#define GWEIGHTS \
const PetscReal gx[4] = {-1.0,  1.0, 1.0, -1.0}, \
                gy[4] = {-1.0, -1.0, 1.0,  1.0};

PetscReal fieldatpt(PetscInt u, PetscInt v, PetscReal xi, PetscReal eta, PetscReal **f) {
  WEIGHTS
  return x[0] * y[0] * f[v][u] + x[1] * y[1] * f[v][u+1] + x[2] * y[2] * f[v+1][u+1] + x[3] * y[3] * f[v+1][u];
}

Grad gradfatpt(PetscInt u, PetscInt v, PetscReal xi, PetscReal eta, PetscReal dx, PetscReal dy, PetscReal **f) {
  Grad gradf;
  WEIGHTS
  GWEIGHTS
  gradf.x =   gx[0] * y[0] * f[v][u]     + gx[1] * y[1] * f[v][u+1]
            + gx[2] * y[2] * f[v+1][u+1] + gx[3] * y[3] * f[v+1][u];
  gradf.y =    x[0] *gy[0] * f[v][u]     +  x[1] *gy[1] * f[v][u+1]
            +  x[2] *gy[2] * f[v+1][u+1] +  x[3] *gy[3] * f[v+1][u];
  gradf.x /= dx;
  gradf.y /= dy;
  return gradf;
}

PetscReal dfieldatpt(PetscInt l, PetscInt u, PetscInt v, PetscReal xi, PetscReal eta) {
  WEIGHTS
  return x[l] * y[l];
}

Grad dgradfatpt(PetscInt l, PetscInt u, PetscInt v, PetscReal xi, PetscReal eta, PetscReal dx, PetscReal dy) {
  Grad dgradfdl;
  WEIGHTS
  GWEIGHTS
  dgradfdl.x = gx[l] *  y[l] / dx;
  dgradfdl.y =  x[l] * gy[l] / dy;
  return dgradfdl;
}
