/* (C) 2016 Ed Bueler */

#include "cmbmodel.h"

PetscErrorCode SetFromOptionsCMBModel(const char *optprefix, CMBModel *cmb) {
  PetscErrorCode ierr;
  cmb->ela   = 2000.0;           // m
  cmb->lapse = 0.001 / secpera;  // s^-1
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,optprefix,"options to cmb model (if used),"");CHKERRQ(ierr);
  FIXME
  ierr = PetscOptionsReal(
      "-D0", "initial (and representative?) value of diffusivity; in m^2/s",
      "continuationscheme.c",cs->D0,&cs->D0,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode GetCMBfromCMBModel(Vec b, Vec H, Vec m, CMBModel *cmb) {
  FIXME
  PetscFunctionReturn(0);
}

