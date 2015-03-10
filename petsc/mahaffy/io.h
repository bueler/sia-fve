/* (C) 2015 Ed Bueler */

#ifndef IO_H_
#define IO_H_

#include <petscvec.h>
#include "mahaffyctx.h"  // for AppCtx

// read special format PETSc binary .dat input file
PetscErrorCode ReadDimensions(AppCtx*);
PetscErrorCode ReadAndReshape2DVec(Vec,PetscViewer,AppCtx*);
PetscErrorCode ReadDataVecs(AppCtx*);

// write PETSc binary .dat files, one per var
PetscErrorCode ViewToBinary(PetscBool,Vec,const char[],const char[]);
PetscErrorCode DumpToFiles(Vec,AppCtx*);

// write a text history file, including command line
PetscErrorCode WriteHistoryFile(const char[],int,char**,AppCtx*);

// use X viewers to show b,m,Hexact
PetscErrorCode ShowFields(AppCtx *user);

#endif // IO_H_
