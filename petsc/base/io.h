/* (C) 2015 Ed Bueler */

/* Methods for reading and writing PETSc binary files, for use with mahaffy.c.
Also utilities for computing stats, viewing at stdout, and viewing with X. */

#ifndef IO_H_
#define IO_H_

#include <petscvec.h>
#include "mahaffyctx.h"  // for AppCtx & DiagnosticScheme

// useful utility; note that VecChop() in PETSc chops any *absolute values*
// below tol
PetscErrorCode VecTrueChop(Vec v, PetscReal tol);

// print according to AppCtx.silent
void myPrintf(const AppCtx *user, const char format[], ...);

// print a report to stdout
PetscErrorCode StdoutReport(Vec,AppCtx *user);

// read a special format PETSc binary .dat input file
PetscErrorCode ReadDimensions(AppCtx *user);
PetscErrorCode ReadAndReshape2DVec(Vec v, PetscViewer viewer, AppCtx *user);
// NOTE: order in which data is read must match order in which nc2petsc.py writes
PetscErrorCode ReadDataVecs(AppCtx *user);

// write a special format PETSc binary .dat output file like the one we read
PetscErrorCode DumpToFile(Vec H, Vec r, AppCtx *user);

// utilities to compute quantities in reports and history files
// volH and volHexact have units m^3; areaH has units m^2
PetscErrorCode GetVolumeArea(Vec H, AppCtx *user, PetscReal *volH, PetscReal *volHexact, PetscReal *areaH);
// enorminf = ||H - Hexact||_infty,  enorm1 = ||H-Hexact||_1;  both have units of m
PetscErrorCode GetErrors(Vec H, AppCtx *user, PetscReal *enorminf, PetscReal *enorm1);
// avD is average and maxD is maximum diffusivity; both have units of m^2 s-1
PetscErrorCode DiffusivityReduce(AppCtx *user, PetscReal *avD, PetscReal *maxD);

// write a text history file, including command line
PetscErrorCode WriteHistoryFile(Vec H, const char name[], int argc, char **argv, AppCtx *user);

// use X viewers to show b,m,Hexact
PetscErrorCode ShowOne(Vec v, PetscInt xdim, PetscInt ydim, const char *title);
PetscErrorCode ShowFields(AppCtx *user);

#endif // IO_H_
