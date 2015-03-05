/* (C) 2015 Ed Bueler */

#ifndef EXACTSIACTX_H_
#define EXACTSIACTX_H_

#include <petscvec.h>
#include "mahaffyctx.h"  // for AppCtx

PetscErrorCode SetToDomeSMB(Vec,PetscBool,const AppCtx*);
PetscErrorCode SetToDomeExactThickness(Vec,const AppCtx*);

/*FIXME
PetscErrorCode SetToJSABed(Vec,const AppCtx*);
PetscErrorCode SetToJSASMB(Vec,const AppCtx*);
PetscErrorCode SetToJSAExactThickness(Vec,const AppCtx*);
*/

#endif // EXACTSIACTX_H_
