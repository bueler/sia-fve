/* (C) 2015 Ed Bueler */

#ifndef BINARYIO_H_
#define BINARYIO_H_

#include <petscvec.h>
#include "mahaffyctx.h"  // for AppCtx

PetscErrorCode ReadDimensions(AppCtx*);
PetscErrorCode ReadAndReshape2DVec(Vec,PetscViewer,AppCtx*);
PetscErrorCode ReadDataVecs(AppCtx*);

PetscErrorCode ViewToBinary(PetscBool,Vec,const char[],const char[]);
PetscErrorCode DumpToFiles(Vec,AppCtx*);

#endif // BINARYIO_H_
