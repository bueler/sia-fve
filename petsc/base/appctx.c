/* (C) 2015 Ed Bueler */

#include <petsc.h>
#include "appctx.h"

PetscErrorCode SetFromOptionsAppCtx(const char *optprefix, AppCtx *user) {
  PetscErrorCode ierr;
  PetscBool      domechosen, dtflg, notryset;
  char           histprefix[512];

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,optprefix,"options to mahaffy","");CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-A", "set value of ice softness A in units Pa-3 s-1",
      NULL,user->A,&user->A,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-averr", "print final average error only, and otherwise run silent",
      NULL,user->averr,&user->averr,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-bedstep", "use bedrock step exact solution by Jarosh, Schoof, Anslow (2013)",
      NULL,user->bedstep,&user->bedstep,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-checkadmissible", "in FormFunctionLocal(), stop if H < 0.0",
      NULL,user->checkadmissible,&user->checkadmissible,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-delta", "dimensionless regularization for slope in SIA formulas",
      NULL,user->delta,&user->delta,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-dome", "use dome exact solution by Bueler (2003) [default]",
      NULL,user->dome,&user->dome,&domechosen);CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-dtjac", "use this time step (years) in Jacobian evaluation, even in steady-state",
      NULL,user->dtjac/user->secpera,&user->dtjac,&dtflg);CHKERRQ(ierr);
  if (dtflg)  user->dtjac *= user->secpera;
  ierr = PetscOptionsReal(
      "-dtrecovery", "use this time step (years) in recovery",
      NULL,user->dtrecovery/user->secpera,&user->dtrecovery,&dtflg);CHKERRQ(ierr);
  if (dtflg)  user->dtrecovery *= user->secpera;
  ierr = PetscOptionsReal(
      "-dt", "use this time step (years) FOR TIME-STEPPING; overrides -mah_dtjac,-mah_dtrecovery",
      NULL,user->dtres/user->secpera,&user->dtres,&dtflg);CHKERRQ(ierr);
  if (dtflg) {
      user->dtres *= user->secpera;
      user->dtjac = user->dtres;
      user->dtrecovery = user->dtres;
  }
  ierr = PetscOptionsString(
      "-dump", "dump fields into PETSc binary files [x,y,b,m,H,Herror].dat with this prefix",
      NULL,user->figsprefix,user->figsprefix,512,&user->dump); CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-initmagic", "constant, in years, used to multiply SMB to get initial iterate for thickness",
      NULL,user->initmagic,&user->initmagic,NULL);CHKERRQ(ierr);
  strcpy(histprefix,"PREFIX/");
  ierr = PetscOptionsString(
      "-history", "write history.txt file with this prefix (also written under -mah_dump)",
      NULL,histprefix,histprefix,512,&user->history); CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-lambda", "amount of upwinding; lambda=0 is none and lambda=1 is full",
      NULL,user->lambda,&user->lambda,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-maxerr", "print final maximum error, and otherwise run silent",
      NULL,user->maxerr,&user->maxerr,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-n", "value of Glen exponent n",
      NULL,user->n,&user->n,NULL);CHKERRQ(ierr);
  if (user->n <= 1.0) {
      SETERRQ1(PETSC_COMM_WORLD,11,"ERROR: n = %f not allowed ... n > 1 is required\n",user->n); }
  ierr = PetscOptionsBool(
      "-nodiag", "do not store, generate, or output diagnostic D and Wmag fields",
      NULL,user->nodiag,&user->nodiag,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-notry", "on SNES diverge, DO NOT try again with recovery",
      NULL,!user->divergetryagain,&user->divergetryagain,&notryset);CHKERRQ(ierr);
  if (notryset) user->divergetryagain = PETSC_FALSE;
  ierr = PetscOptionsString(
      "-read", "read grid and data from special-format PETSc binary file; see README.md",
      NULL,user->readname,user->readname,512,&user->read); CHKERRQ(ierr);
  ierr = PetscOptionsString(
      "-readinitial", "read grid initial H from this special-format PETSc binary file",
      NULL,user->readinitialname,user->readinitialname,512,&user->readinitial); CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-showdata", "use PETSc X viewers to show b, m, and exact/observed H",
      NULL,user->showdata,&user->showdata,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-silent", "run silent (print nothing)",
      NULL,user->silent,&user->silent,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt(
      "-steps", "number of Backward Euler time-steps to take; default=1 is for steady-state",
      NULL,user->numsteps,&user->numsteps,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-swapxy", "swap coordinates x and y when building bedrock step exact solution",
      NULL,user->swapxy,&user->swapxy,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-true", "use true Mahaffy method, not default M*",
      NULL,user->mtrue,&user->mtrue,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  // enforce consistency of cases
  if (user->numsteps > 1) {
      if ((user->dtres <= 0.0) || (user->dtjac <= 0.0) || (user->dtres != user->dtjac)) {
        SETERRQ(PETSC_COMM_WORLD,4,
           "OPTION CONFLICT: Backward Euler time-steps requested but dtres and dtjac are inconsistent; -mah_dt not set?\n");
      }
  }
  if ((user->averr) || (user->maxerr))
      user->silent = PETSC_TRUE;
  if (user->read) {
      if (domechosen) {
        SETERRQ(PETSC_COMM_WORLD,1,"OPTION CONFLICT: both -mah_dome and -mah_read not allowed\n");
      } else
        user->dome = PETSC_FALSE;    // here user has chosen -mah_read and not -mah_dome
  }
  if (user->bedstep) {
      if (domechosen) {
        SETERRQ(PETSC_COMM_WORLD,2,"ERROR option conflict: both -mah_dome and -mah_bedstep not allowed\n");
      } else
        user->dome = PETSC_FALSE;    // here user has chosen -mah_bedstep and not -mah_dome
  }
  if ((user->read) && (user->bedstep)) {
      SETERRQ(PETSC_COMM_WORLD,3,"ERROR option conflict: both -mah_bedstep and -mah_read not allowed\n");
  }
  if (user->mtrue)
      user->lambda = 0.0; // force no upwinding
  if (user->dump)
      user->history = PETSC_TRUE; // history is part of dump
  else if (user->history)
      strcpy(user->figsprefix,histprefix); // store path to history file in figsprefix

  PetscFunctionReturn(0);
}

