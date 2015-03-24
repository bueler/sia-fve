/* (C) 2015 Ed Bueler */

/* Formulas for two exact solutions, one from Bueler2003 and one from JaroschSchoofAnslow2013. */

#ifndef EXACTSIACTX_H_
#define EXACTSIACTX_H_

#include <petscvec.h>
#include "mahaffyctx.h"  // for AppCtx

PetscErrorCode DomeCMB(Vec,const AppCtx*);
PetscErrorCode DomeExactThickness(Vec,const AppCtx*);

PetscErrorCode BedStepBed(Vec,const AppCtx*);
PetscErrorCode BedStepCMB(Vec,const AppCtx*);
PetscErrorCode BedStepExactThickness(Vec,const AppCtx*);

#endif // EXACTSIACTX_H_
