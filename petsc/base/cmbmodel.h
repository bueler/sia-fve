/* (C) 2016 Ed Bueler */

#ifndef CMBMODEL_H_
#define CMBMODEL_H_

#include <petsc.h>

// model is:
//     m(x,y) = lapse * (s(x,y) - ela)
// where m(x,y) in m s^-1 is the climatic mass balance and where
//     s(x,y) = b(x,y) + H(x,y)
// is the surface elevation in m

typedef struct {
  PetscReal ela,   // elevation in meters of equilibrium-line-altitude
            lapse; // lapse rate in s^-1
} CMBModel;

PetscErrorCode SetFromOptionsCMBModel(CMBModel *cmb, const char *optprefix, PetscReal secpera);

PetscErrorCode M_CMBModel(CMBModel *cmb, PetscReal b, PetscReal H, PetscReal *m);

PetscErrorCode dMdH_CMBModel(CMBModel *cmb, PetscReal *dmdH);

#endif // CMBMODEL_H_

