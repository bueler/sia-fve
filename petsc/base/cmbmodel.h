/* (C) 2016 Ed Bueler */

#ifndef CMBMODEL_H_
#define CMBMODEL_H_

#include <petsc.h>

// model is:
//     M = cmbgrad * (s - ela)
// where  M  in m s^-1 is the climatic mass balance and where
// s  is the surface elevation in m

typedef struct {
  PetscReal ela,   // equilibrium line altitude (m)
            zgrad; // vertical derivative (gradient) of CMB (s^-1)
} CMBModel;

PetscErrorCode SetFromOptionsCMBModel(CMBModel *cmb, const char *optprefix, PetscReal secpera);

PetscErrorCode M_CMBModel(CMBModel *cmb, PetscReal s, PetscReal *M);

PetscErrorCode dMds_CMBModel(CMBModel *cmb, PetscReal s, PetscReal *dMds);

#endif // CMBMODEL_H_

