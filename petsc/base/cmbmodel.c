/* (C) 2016 Ed Bueler */

#include "cmbmodel.h"

PetscErrorCode SetFromOptionsCMBModel(CMBModel *cmb, const char *optprefix, PetscReal secpera) {
  PetscErrorCode ierr;
  cmb->ela   = 2000.0; // m
  cmb->zgrad = 0.001;  // a^-1
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,optprefix,
            "options to climatic mass balance (CMB) model, if used","");CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-ela", "equilibrium line altitude, in m",
      "cmbmodel.c",cmb->ela,&cmb->ela,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-zgrad", "vertical derivative (gradient) of CMB, in a^-1",
      "cmbmodel.c",cmb->zgrad,&cmb->zgrad,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  cmb->zgrad /= secpera;
  PetscFunctionReturn(0);
}


PetscErrorCode M_CMBModel(CMBModel *cmb, PetscReal s, PetscReal *M) {
  *M = cmb->zgrad * (s - cmb->ela);
  PetscFunctionReturn(0);
}


PetscErrorCode dMdH_CMBModel(CMBModel *cmb, PetscReal s, PetscReal *dMds) {
  *dMds = cmb->zgrad;
  PetscFunctionReturn(0);
}

