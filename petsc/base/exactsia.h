/* (C) 2015 Ed Bueler */

/* Formulas for two exact solutions, one from Bueler2003 and one from JaroschSchoofAnslow2013. */

#ifndef EXACTSIA_H_
#define EXACTSIA_H_

#include <petscvec.h>
#include "appctx.h"

PetscErrorCode DomeCMB(Vec,const AppCtx*);
PetscErrorCode DomeExactThickness(Vec,const AppCtx*);

PetscErrorCode BedStepBed(Vec,const AppCtx*);
PetscErrorCode BedStepCMB(Vec,const AppCtx*);
PetscErrorCode BedStepExactThickness(Vec,const AppCtx*);

#endif // EXACTSIA_H_
