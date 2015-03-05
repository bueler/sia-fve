/* (C) 2015 Ed Bueler */

#ifndef EXACTSIACTX_H_
#define EXACTSIACTX_H_

#include <petscvec.h>
#include "mahaffyctx.h"  // for AppCtx

PetscErrorCode SetToDomeCMB(Vec,PetscBool,const AppCtx*);
PetscErrorCode SetToDomeExactThickness(Vec,const AppCtx*);

PetscErrorCode SetToBedStepBed(Vec,const AppCtx*);
PetscErrorCode SetToBedStepCMB(Vec,const AppCtx*);
PetscErrorCode SetToBedStepExactThickness(Vec,const AppCtx*);

#endif // EXACTSIACTX_H_
