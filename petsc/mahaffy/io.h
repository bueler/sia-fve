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

// utilities to compute quantities in reports and history files
PetscErrorCode GetMaxDiffusivity(AppCtx*,PetscReal*);
PetscErrorCode GetVolumes(Vec,AppCtx*,PetscReal*,PetscReal*);
PetscErrorCode GetErrors(Vec,AppCtx*,PetscReal*,PetscReal*);

// print a report to stdout
PetscErrorCode StdoutReport(Vec,DMDALocalInfo*,AppCtx*);

// write a text history file, including command line
PetscErrorCode WriteHistoryFile(Vec,const char[],int,char**,DMDALocalInfo*,AppCtx*);

// use X viewers to show b,m,Hexact
PetscErrorCode ShowFields(AppCtx*);

#endif // IO_H_
