/* (C) 2016 Ed Bueler */

#include "cmbmodel.h"

PetscErrorCode SetFromOptionsCMBModel(CMBModel *cmb, const char *optprefix, PetscReal secpera) {
  PetscErrorCode ierr;
  cmb->ela   = 2000.0; // m
  cmb->lapse = 0.001;  // a^-1
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,optprefix,
            "options to climatic mass balance (CMB) model, if used","");CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-ela", "equilibrium line altitude in m",
      "cmbmodel.c",cmb->ela,&cmb->ela,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-lapse", "lapse rate, the rate of CMB increase per meter of elevation gain, in a^-1",
      "cmbmodel.c",cmb->lapse,&cmb->lapse,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  cmb->lapse /= secpera;
  PetscFunctionReturn(0);
}


PetscErrorCode M_CMBModel(CMBModel *cmb, PetscReal b, PetscReal H, PetscReal *m) {
  *m = cmb->lapse * (b + H - cmb->ela);
  PetscFunctionReturn(0);
}


PetscErrorCode dMdH_CMBModel(CMBModel *cmb, PetscReal *dmdH) {
  *dmdH = cmb->lapse;
  PetscFunctionReturn(0);
}

