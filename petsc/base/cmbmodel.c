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


PetscErrorCode M_CMBModel(CMBModel *cmb, DM da, Vec b, Vec H, Vec m) {
  PetscErrorCode ierr;
  DMDALocalInfo info;
  PetscReal     **am, **ab, **aH;
  PetscInt      j, k;

  PetscFunctionBeginUser;
  ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, m, &am);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, b, &ab);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, H, &aH);CHKERRQ(ierr);
  for (k=info.ys; k<info.ys+info.ym; k++) {
      for (j=info.xs; j<info.xs+info.xm; j++) {
          am[k][j] = cmb->lapse * (ab[k][j] + aH[k][j] - cmb->ela);
      }
  }
  ierr = DMDAVecRestoreArray(da, m, &am);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, b, &ab);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, H, &aH);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode dMdH_CMBModel(CMBModel *cmb, PetscReal *dmdH) {
  *dmdH = cmb->lapse;
  PetscFunctionReturn(0);
}

