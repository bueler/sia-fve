/* (C) 2015 Ed Bueler */

#include "q1op.h"

PetscErrorCode fieldatpt(PetscInt j, PetscInt k,         // (j,k) is the element (by lower-left corner)
                         PetscReal locx, PetscReal locy, // = (x,y) coords in element
                         PetscReal **f,                  // f[k][j] are node values
                         const AppCtx *user, PetscReal *f_atpt) {
  const PetscReal dx = user->dx,  dy = dx,
                  x[4]  = {1.0 - locx / dx, locx / dx,       1.0 - locx / dx, locx / dx},
                  y[4]  = {1.0 - locy / dy, 1.0 - locy / dy, locy / dy,       locy / dy};
  if (f == NULL)  SETERRQ(PETSC_COMM_WORLD,1,"ERROR: illegal NULL ptr in fieldatpt() ...\n");
  *f_atpt = x[0] * y[0] * f[k][j] + x[1] * y[1] * f[k][j+1] + x[2] * y[2] * f[k+1][j] + x[3] * y[3] * f[k+1][j+1];
  PetscFunctionReturn(0);
}

PetscErrorCode gradatpt(PetscInt j, PetscInt k,         // (j,k) is the element (by lower-left corner)
                        PetscReal locx, PetscReal locy, // = (x,y) coords in element
                        PetscReal **f,                  // f[k][j] are node values
                        const AppCtx *user, Grad *gradf_atpt) {
  const PetscReal dx = user->dx, dy = dx,
                  x[4]  = {1.0 - locx / dx, locx / dx,       1.0 - locx / dx, locx / dx},
                  gx[4] = {- 1.0 / dx,      1.0 / dx,        - 1.0 / dx,      1.0 / dx},
                  y[4]  = {1.0 - locy / dy, 1.0 - locy / dy, locy / dy,       locy / dy},
                  gy[4] = {- 1.0 / dy,      - 1.0 / dy,      1.0 / dy,        1.0 / dy};
  if (f == NULL)  SETERRQ(PETSC_COMM_WORLD,1,"ERROR: illegal NULL ptr in gradatpt() ...\n");
  gradf_atpt->x =   gx[0] * y[0] * f[k][j]   + gx[1] * y[1] * f[k][j+1]
                  + gx[2] * y[2] * f[k+1][j] + gx[3] * y[3] * f[k+1][j+1];
  gradf_atpt->y =    x[0] *gy[0] * f[k][j]   +  x[1] *gy[1] * f[k][j+1]
                  +  x[2] *gy[2] * f[k+1][j] +  x[3] *gy[3] * f[k+1][j+1];
  PetscFunctionReturn(0);
}
