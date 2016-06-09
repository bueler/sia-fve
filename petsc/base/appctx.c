/* (C) 2015 Ed Bueler */

#include <petsc.h>
#include "appctx.h"

// default settings of parameters
PetscErrorCode initialize(AppCtx *user) {
  user->n      = 3.0;
  user->g      = 9.81;       // m/s^2
  user->rho    = 910.0;      // kg/m^3
  user->secpera= 31556926.0;
  user->A      = 1.0e-16/user->secpera; // = 3.17e-24  1/(Pa^3 s); EISMINT I value

  user->initmagic = 1000.0;  // a
  user->delta  = 1.0e-4;

  user->lambda = 0.25;  // amount of upwinding; some trial-and-error with bedstep soln; 0.1 gives some Newton convergence difficulties on refined grid (=125m); earlier M* used 0.5

  user->T          = -1.0;
  user->dtres      = 0.0;
  user->dtjac      = 0.0;
  user->dtrecovery = 1.0 * user->secpera;  // default 1 year time step for Backward Euler
  user->goodm      = -1;
  user->recoverycount = 0;

  user->dumpdt     = -user->secpera;

  user->mtrue      = PETSC_FALSE;

  user->dome       = PETSC_TRUE;  // defaults to this case
  user->bedstep    = PETSC_FALSE;
  user->swapxy     = PETSC_FALSE;
  user->divergetryagain = PETSC_TRUE;
  user->checkadmissible = PETSC_FALSE;

  user->read       = PETSC_FALSE;
  user->readinitial= PETSC_FALSE;
  user->showdata   = PETSC_FALSE;
  user->history    = PETSC_FALSE;
  user->nodiag     = PETSC_FALSE;
  user->dump       = PETSC_FALSE;
  user->silent     = PETSC_FALSE;
  user->averr      = PETSC_FALSE;
  user->maxerr     = PETSC_FALSE;
  user->cmbmodel   = PETSC_FALSE;

  strcpy(user->figsprefix,"PREFIX/");  // dummies improve "mahaffy -help" appearance
  strcpy(user->readname,"FILENAME");
  strcpy(user->readinitialname,"FILENAME");

  user->cs = NULL;
  user->cmb = NULL;

  PetscFunctionReturn(0);
}


PetscErrorCode SetFromOptionsAppCtx(const char *optprefix, AppCtx *user) {
  PetscErrorCode ierr;
  PetscBool      domechosen, dtflg, notryset, Tset;
  char           histprefix[512], initHname[512], initsname[512];

  ierr = initialize(user); CHKERRQ(ierr);
  strcpy(initHname,"FILENAME");
  strcpy(initsname,"FILENAME");
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,optprefix,"options to mahaffy","");CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-A", "set value of ice softness A in units Pa-3 s-1",
      "mahaffy.c",user->A,&user->A,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-averr", "print final average error only, and otherwise run silent",
      "mahaffy.c",user->averr,&user->averr,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-bedstep", "use bedrock step exact solution by Jarosh, Schoof, Anslow (2013)",
      "mahaffy.c",user->bedstep,&user->bedstep,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-checkadmissible", "in FormFunctionLocal(), stop if H < 0.0",
      "mahaffy.c",user->checkadmissible,&user->checkadmissible,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-cmbmodel", "use ela + lapse-rate model for computing climatic mass balance (CMB)",
      "mahaffy.c",user->cmbmodel,&user->cmbmodel,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-delta", "dimensionless regularization for slope in SIA formulas",
      "mahaffy.c",user->delta,&user->delta,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-dome", "use dome exact solution by Bueler (2003) [default]",
      "mahaffy.c",user->dome,&user->dome,&domechosen);CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-dtjac", "in steady-state, use this time step (years) in Jacobian evaluation",
      "mahaffy.c",user->dtjac/user->secpera,&user->dtjac,&dtflg);CHKERRQ(ierr);
  if (dtflg)  user->dtjac *= user->secpera;
  ierr = PetscOptionsReal(
      "-dtrecovery", "in steady-state, use this time step (years) in recovery",
      "mahaffy.c",user->dtrecovery/user->secpera,&user->dtrecovery,&dtflg);CHKERRQ(ierr);
  if (dtflg)  user->dtrecovery *= user->secpera;
  ierr = PetscOptionsReal(
      "-dt", "use this time step (years); overrides -mah_dtjac,-mah_dtrecovery",
      "mahaffy.c",user->dtres/user->secpera,&user->dtres,&dtflg);CHKERRQ(ierr);
  if (dtflg) {
      user->dtres *= user->secpera;
      user->dtjac = user->dtres;
      user->dtrecovery = user->dtres;
  }
  ierr = PetscOptionsString(
      "-dump", "dump final fields into PETSc binary files [x,y,b,m,H,Herror].dat with this prefix",
      "mahaffy.c",user->figsprefix,user->figsprefix,512,&user->dump); CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-dumpdt", "period between dumping states when time-stepping, in years; setting positive activates",
      "mahaffy.c",user->dumpdt/user->secpera,&user->dumpdt,&dtflg);CHKERRQ(ierr);
  if (dtflg)  user->dumpdt *= user->secpera;
  ierr = PetscOptionsReal(
      "-initmagic", "constant, in years, used to multiply CMB to get initial iterate for thickness",
      "mahaffy.c",user->initmagic,&user->initmagic,NULL);CHKERRQ(ierr);
  strcpy(histprefix,"PREFIX/");
  ierr = PetscOptionsString(
      "-history", "write history.txt file with this prefix (also written under -mah_dump)",
      "mahaffy.c",histprefix,histprefix,512,&user->history); CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-lambda", "amount of upwinding; lambda=0 is none and lambda=1 is full",
      "mahaffy.c",user->lambda,&user->lambda,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-maxerr", "print final maximum error, and otherwise run silent",
      "mahaffy.c",user->maxerr,&user->maxerr,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-n", "value of Glen exponent n",
      "mahaffy.c",user->n,&user->n,NULL);CHKERRQ(ierr);
  if (user->n <= 1.0) {
      SETERRQ1(PETSC_COMM_WORLD,11,"ERROR: n = %f not allowed ... n > 1 is required\n",user->n); }
  ierr = PetscOptionsBool(
      "-nodiag", "do not store, generate, or output diagnostic D and Wmag fields",
      "mahaffy.c",user->nodiag,&user->nodiag,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-notry", "in steady-state, on SNES diverge, DO NOT try again with recovery",
      "mahaffy.c",!user->divergetryagain,&user->divergetryagain,&notryset);CHKERRQ(ierr);
  if (notryset) user->divergetryagain = PETSC_FALSE;
  ierr = PetscOptionsString(
      "-read", "read grid and data from special-format PETSc binary file; see README.md",
      "mahaffy.c",user->readname,user->readname,512,&user->read); CHKERRQ(ierr);
  ierr = PetscOptionsString(
      "-readinitial", "read initial thickness H from this special-format PETSc binary file",
      "mahaffy.c",initHname,initHname,512,&user->readinitial); CHKERRQ(ierr);
  ierr = PetscOptionsString(
      "-readinitialsurface", "generate initial H by reading surface s from this file",
      "mahaffy.c",initsname,initsname,512,&user->readinitialsurface); CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-showdata", "use PETSc X viewers to show b, m, and exact/observed H",
      "mahaffy.c",user->showdata,&user->showdata,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-silent", "run silent (print nothing)",
      "mahaffy.c",user->silent,&user->silent,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-swapxy", "swap coordinates x and y when building bedrock step exact solution",
      "mahaffy.c",user->swapxy,&user->swapxy,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-T", "use Backward Euler time-stepping to reach this time (years)",
      "mahaffy.c",user->T/user->secpera,&user->T,&Tset);CHKERRQ(ierr);
  if (Tset)  user->T *= user->secpera;
  ierr = PetscOptionsBool(
      "-true", "use true Mahaffy method, not default M*",
      "mahaffy.c",user->mtrue,&user->mtrue,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  // enforce consistency of cases
  if (Tset) {
      if ((user->dtres <= 0.0) || (user->dtjac <= 0.0) || (user->dtres != user->dtjac)) {
        SETERRQ(PETSC_COMM_WORLD,4,
           "OPTION CONFLICT: Backward Euler time-steppping requested (-mah_T set) but dtres and dtjac are inconsistent; -mah_dt not set?\n");
      }
      if (user->T <= 0) {
          SETERRQ(PETSC_COMM_WORLD,6,"OPTION CONFLICT: -mah_T value must be positive\n");
      }
      if (user->T < user->dtres) {
          SETERRQ(PETSC_COMM_WORLD,7,"OPTION CONFLICT: -mah_T value must be at least as large as -mah_dt value\n");
      }
  } else if (user->dtres > 0.0) {
      SETERRQ(PETSC_COMM_WORLD,5,"OPTION CONFLICT: dtres > 0  but no time-stepping not requested; set -mah_T?\n");
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
  if (user->readinitial) {
      if (user->readinitialsurface) {
          SETERRQ(PETSC_COMM_WORLD,8,
              "ERROR option conflict: both -mah_readinitial and -mah_readinitialsurface not allowed\n");
      }
      strcpy(user->readinitialname,initHname);
  } else if (user->readinitialsurface)
      strcpy(user->readinitialname,initsname);
  if (user->mtrue)
      user->lambda = 0.0; // force no upwinding if user wants true Mahaffy
  if (user->dump)
      user->history = PETSC_TRUE; // history is part of dump
  else if (user->history)
      strcpy(user->figsprefix,histprefix); // store path to history file in figsprefix

  // derived constant computed after n,A get set
  user->Gamma = 2.0 * PetscPowReal(user->rho*user->g,user->n) * user->A / (user->n+2.0);

  PetscFunctionReturn(0);
}

