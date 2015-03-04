/* (C) 2015 Ed Bueler */

#ifndef EXACTSIACTX_H_   /* Include guard */
#define EXACTSIACTX_H_

#include <petscdmda.h>
#include <petscvec.h>

#include "mahaffyctx.h"

PetscErrorCode SetToDomeSMB(Vec,PetscBool,const AppCtx*);
PetscErrorCode SetToDomeExactThickness(Vec,const AppCtx*);

/*FIXME
PetscErrorCode SetToJSABed(Vec,const AppCtx*);
PetscErrorCode SetToJSASMB(Vec,const AppCtx*);
PetscErrorCode SetToJSAExactThickness(Vec,const AppCtx*);
*/

#endif // EXACTSIACTX_H_
