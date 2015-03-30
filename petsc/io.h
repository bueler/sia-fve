/* (C) 2015 Ed Bueler */

/* Methods for reading and writing PETSc binary files, for use with mahaffy.c.
Also utilities for computing stats, viewing at stdout, and viewing with X. */

#ifndef IO_H_
#define IO_H_

#include <petscvec.h>
#include "mahaffyctx.h"  // for AppCtx

// print according to AppCtx.silent
void myPrintf(const AppCtx*, const char[], ...);

// print a report to stdout
PetscErrorCode StdoutReport(Vec,AppCtx*);

// read special format PETSc binary .dat input file
PetscErrorCode ReadDimensions(AppCtx*);
PetscErrorCode ReadAndReshape2DVec(Vec,PetscViewer,AppCtx*);
PetscErrorCode ReadDataVecs(AppCtx*);

// write PETSc binary .dat files, one per var
PetscErrorCode ViewToBinary(PetscBool,Vec,const char[],const char[]);
PetscErrorCode DumpToFiles(Vec,AppCtx*);

// utilities to compute quantities in reports and history files
PetscErrorCode GetVolumes(Vec,AppCtx*,PetscReal*,PetscReal*);
PetscErrorCode GetErrors(Vec,AppCtx*,PetscReal*,PetscReal*);

// write a text history file, including command line
PetscErrorCode WriteHistoryFile(Vec,const char[],int,char**,AppCtx*);

// use X viewers to show b,m,Hexact
PetscErrorCode ShowFields(AppCtx*);

#endif // IO_H_
