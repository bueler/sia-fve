/* (C) 2015 Ed Bueler */

#ifndef MAHIO_H_
#define MAHIO_H_

#include <petscvec.h>
#include "mahaffyctx.h"  // for AppCtx

typedef struct {
  Vec  topgread, cmbread, thkobsread;  // only valid if user.read is TRUE
} SerialReadVecs;

PetscErrorCode ViewToBinary(PetscBool,Vec,const char[],const char[]);
PetscErrorCode DumpToFiles(Vec,AppCtx*);
PetscErrorCode CreateReadVecs(SerialReadVecs*,AppCtx*);
PetscErrorCode ReshapeAndDestroyReadVecs(SerialReadVecs*,AppCtx*);

#endif // MAHIO_H_
