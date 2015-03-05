/* (C) 2015 Ed Bueler */

#ifndef BINARYIO_H_
#define BINARYIO_H_

#include <petscvec.h>
#include "mahaffyctx.h"  // for AppCtx

// read special format PETSc binary .dat input file
PetscErrorCode ReadDimensions(AppCtx*);
PetscErrorCode ReadAndReshape2DVec(Vec,PetscViewer,AppCtx*);
PetscErrorCode ReadDataVecs(AppCtx*);

// write PETSc binary .day files, one per var
PetscErrorCode ViewToBinary(PetscBool,Vec,const char[],const char[]);
PetscErrorCode DumpToFiles(Vec,AppCtx*);

// use X viewers to show b,m,Hexact
PetscErrorCode ShowFields(AppCtx *user);

#endif // BINARYIO_H_
