/* (C) 2015 Ed Bueler */

static const char help[] =
"Solves steady ice sheet problem in 2d:\n"
"    div (q^x,q^y) = m,\n"
"    (q^x,q^y) = - Gamma H^{n+2} |grad s|^{n-1} grad s,\n"
"where  H(x,y)  is the ice thickness,  b(x,y)  is the bed elevation,\n"
"and  s(x,y) = H(x,y) + b(x,y)  is surface elevation.\n"
"Note  n > 1  and  Gamma = 2 A (rho g)^n / (n+2).\n"
"Domain is  -Lx < x < Lx,  -Ly < y < Ly,  with periodic boundary conditions.\n\n"
"Computed by Q1 FVE method with either analytical or FD evaluation of Jacobian.\n"
"Compares M* improved scheme (default), different amounts of upwinding, and\n"
"true Mahaffy schemes.  Uses SNESVI with constraint  H(x,y) >= 0.\n\n"
"Three problem cases:\n"
"  (1) flat bed case where analytical solution is known\n"
"  (2) bedrock-step case where analytical solution is known\n"
"  (3) reads input data b and m from PETSc binary file; see grn/README.md.\n\n";

/* See README.md first.  Here are additional technical notes:

Solution options:
   ./mahaffy -cs_end 4         # don't go all the way on continuation scheme
   ./mahaffy -cs_D0 1.0        # change constant diffusivity used in continuation scheme to other value than
                               # default = 10.0;  big may be good for convergence, esp. w upwinding
   ./mahaffy -mah_bedstep -mah_lambda 0.0  # NO upwinding on  grad b  part of flux
   ./mahaffy -mah_bedstep -mah_lambda 1.0  # FULL upwinding
   ./mahaffy -mah_notry        # on SNES diverge, do not try again in recovery mode

Jacobian by FD coloring:
   ./mahaffy -snes_fd_color    # with M*; this requires 3|mx and 3|my !

True Mahaffy, i.e. M* sans quadrature & upwind improvements, *only* runs with
Jacobian by FD coloring:
   ./mahaffy -mah_true -snes_fd_color -da_grid_x 15 -da_grid_y 15   # this requires 5|mx and 5|my !

PETSc solver variations:
   ./mahaffy -snes_type vinewtonssls  # vinewtonrsls is default

Fully converges for these levels:
   for LEV in 0 1 2 3; do  mpiexec -n 6 ./mahaffy -da_refine $LEV -snes_type vinewtonssls -snes_max_it 200 -pc_type asm -sub_pc_type lu -mah_notry; done

Divergence:
   ./mahaffy -da_refine 2                   # DIVERGED_LINE_SEARCH at 10
   ./mahaffy -da_refine 2 -snes_type vinewtonssls # DIVERGED       at 12
   ./mahaffy -mah_bedstep -mah_notry        # DIVERGED_LU_ZERO_PIVOT at 12

Do "manual" time steps after divergence, to approach steady state:
   mkdir teststeps/
   mpiexec -n 6 ./mahaffy -da_refine 2 -mah_dump teststeps/ -mah_notry
   # above diverges at 10 with        errors:  max =  489.93 m,  av =   24.21 m,  voldiff% =  2.10
   mpiexec -n 6 ./mahaffy -mah_read teststeps/unnamed.dat -mah_readinitial -cs_start 8 -mah_steps 100 -mah_dtBE 100.0 -mah_notry
   # above completes 10000 a run with errors:  max =  412.68 m,  av =    6.52 m,  voldiff% =  0.34
*/

#include <math.h>
#include <petscdmda.h>
#include <petscsnes.h>

#include "base/mahaffyctx.h"
#include "base/continuationscheme.h"
#include "base/exactsia.h"
#include "base/io.h"
#include "base/solver.h"

extern PetscErrorCode FormBounds(SNES,Vec,Vec);
extern PetscErrorCode ProcessOptions(AppCtx*);
extern PetscErrorCode ChopScaleSMBforInitial(Vec,AppCtx*);


int main(int argc,char **argv) {
  PetscErrorCode      ierr;
  SNES                snes;
  Vec                 H, Htry;
  AppCtx              user;
  DMDALocalInfo       info;
  PetscReal           dx,dy;
  PetscInt            l, m;

  PetscInitialize(&argc,&argv,(char*)0,help);

  // default settings of parameters
  user.n      = 3.0;
  user.g      = 9.81;       // m/s^2
  user.rho    = 910.0;      // kg/m^3
  user.secpera= 31556926.0;
  user.A      = 1.0e-16/user.secpera; // = 3.17e-24  1/(Pa^3 s); EISMINT I value

  user.initmagic = 1000.0;  // a
  user.delta  = 1.0e-4;

  user.lambda = 0.25;  // amount of upwinding; some trial-and-error with bedstep soln; 0.1 gives some Newton convergence difficulties on refined grid (=125m); earlier M* used 0.5

  user.doBEsteps  = PETSC_FALSE;
  user.numBEsteps = 1;                   // recovery mode does one step
  user.dtBE       = 1.0 * user.secpera;  // default 1 year time step for Backward Euler

  user.recoverycount = 0;

  user.mtrue      = PETSC_FALSE;

  user.dome       = PETSC_TRUE;  // defaults to this case
  user.bedstep    = PETSC_FALSE;
  user.swapxy     = PETSC_FALSE;
  user.divergetryagain = PETSC_TRUE;
  user.checkadmissible = PETSC_FALSE;

  user.read       = PETSC_FALSE;
  user.readinitial= PETSC_FALSE;
  user.showdata   = PETSC_FALSE;
  user.history    = PETSC_FALSE;
  user.nodiag     = PETSC_FALSE;
  user.dump       = PETSC_FALSE;

  user.silent     = PETSC_FALSE;
  user.averr      = PETSC_FALSE;
  user.maxerr     = PETSC_FALSE;

  ierr = InitializeCS(&(user.cs)); CHKERRQ(ierr);

  strcpy(user.figsprefix,"PREFIX/");  // dummies improve "mahaffy -help" appearance
  strcpy(user.readname,"FILENAME");

  ierr = ProcessOptions(&user); CHKERRQ(ierr);

  // derived constant computed after n,A get set
  user.Gamma = 2.0 * PetscPowReal(user.rho*user.g,user.n) * user.A / (user.n+2.0);

  if (user.read) {
      myPrintf(&user,
          "reading dimensions and grid spacing from %s ...\n",user.readname);
      ierr = ReadDimensions(&user); CHKERRQ(ierr);  // sets user.[Nx,Ny,Lx,Ly,dx]
  } else {
      if (user.dome) {
          user.Nx = -18;        // so DMDACreate2d() defaults to 18x18
          user.Ny = -18;
          user.Lx = 900.0e3;    // m
          user.Ly = 900.0e3;    // m
      } else {
          user.Nx = -30;
          user.Ny = -30;
          user.Lx = 30.0e3;    // m
          user.Ly = 30.0e3;    // m
      }
  }

  // this DMDA is used for scalar fields on nodes; cell-centered grid
  ierr = DMDACreate2d(PETSC_COMM_WORLD,
                      DM_BOUNDARY_PERIODIC,DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX,
                      user.Nx,user.Ny,PETSC_DECIDE,PETSC_DECIDE,  // default grid if Nx<0, Ny<0
                      1, (user.mtrue==PETSC_TRUE) ? 2 : 1,        // dof=1, stencilwidth=1
                      NULL,NULL,&user.da);
  ierr = DMSetApplicationContext(user.da, &user);CHKERRQ(ierr);

  // we now know grid spacing
  ierr = DMDAGetLocalInfo(user.da,&info); CHKERRQ(ierr);
  if (user.read == PETSC_FALSE) {
      user.dx = 2.0 * user.Lx / (PetscReal)(info.mx);
  }
  dx = user.dx; // for now,
  dy = user.dx; // square elements
  ierr = DMDASetUniformCoordinates(user.da, -user.Lx+dx/2, user.Lx+dx/2, -user.Ly+dy/2, user.Ly+dy/2,
                                   0.0,1.0); CHKERRQ(ierr);
  myPrintf(&user,"solving on [-Lx,Lx]x[-Ly,Ly] with  Lx=%.3f km  and  Ly=%.3f km\n"
                 "grid of  %d x %d  points with spacing  dx = %.6f km ...\n",
           user.Lx/1000.0,user.Ly/1000.0,info.mx,info.my,dx/1000.0);

  // this DMDA is used for evaluating fluxes at 4 quadrature points on each element
  ierr = DMDACreate2d(PETSC_COMM_WORLD,
                      DM_BOUNDARY_PERIODIC,DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX,
                      info.mx,info.my,PETSC_DECIDE,PETSC_DECIDE,
                      4, (user.mtrue==PETSC_TRUE) ? 2 : 1,  // dof=4,  stencilwidth=1 or 2
                      NULL,NULL,&user.quadda);
  ierr = DMSetApplicationContext(user.quadda, &user);CHKERRQ(ierr);

  // this DMDA is used for evaluating DfluxDl at 4 quadrature points on each
  // elements but with respect to 4 nodal values
  ierr = DMDACreate2d(PETSC_COMM_WORLD,
                      DM_BOUNDARY_PERIODIC,DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX,
                      info.mx,info.my,PETSC_DECIDE,PETSC_DECIDE,
                      16, 1,  // dof=16,  stencilwidth=1 ALWAYS
                              // SETERRQ() in Jacobian protects from use in mtrue=TRUE case
                      NULL,NULL,&user.sixteenda);
  ierr = DMSetApplicationContext(user.sixteenda, &user);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(user.da,&H);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)H,"thickness solution H"); CHKERRQ(ierr);

  ierr = VecDuplicate(H,&Htry); CHKERRQ(ierr);
  ierr = VecDuplicate(H,&user.Hinitial); CHKERRQ(ierr);

  ierr = VecDuplicate(H,&user.b); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)(user.b),"bed elevation b"); CHKERRQ(ierr);
  ierr = VecDuplicate(H,&user.m); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)(user.m),"surface mass balance m"); CHKERRQ(ierr);
  ierr = VecDuplicate(H,&user.Hexact); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)(user.Hexact),"exact/observed thickness H"); CHKERRQ(ierr);

  ierr = VecDuplicate(H,&user.ds.Dnodemax); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)(user.ds.Dnodemax),"maximum diffusivity D at node"); CHKERRQ(ierr);
  ierr = VecDuplicate(H,&user.ds.Wmagnodemax); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)(user.ds.Wmagnodemax),"maximum pseudo-velocity W magnitude at node"); CHKERRQ(ierr);

  // fill user.[b,m,Hexact] according to 3 choices: data, dome exact soln, JSA exact soln
  if (user.read) {
      myPrintf(&user,"reading b, m, Hexact (or Hobserved), Hinitial from %s ...\n", user.readname);
      ierr = ReadDataVecs(&user); CHKERRQ(ierr);
      if (!user.readinitial) {
          myPrintf(&user,"  ignoring read Hinitial ...\n");
          ierr = ChopScaleSMBforInitial(user.Hinitial,&user); CHKERRQ(ierr);
      }
  } else {
      if (user.dome) {
          myPrintf(&user,"generating b, m, Hexact from dome formulas in Bueler (2003) ...\n");
          ierr = VecSet(user.b,0.0); CHKERRQ(ierr);
          ierr = DomeCMB(user.m,&user); CHKERRQ(ierr);
          ierr = DomeExactThickness(user.Hexact,&user); CHKERRQ(ierr);
      } else if (user.bedstep) {
          myPrintf(&user,"generating b, m, Hexact from bedrock step formulas in Jarosch et al (2013) ...\n");
          ierr = BedStepBed(user.b,&user); CHKERRQ(ierr);
          ierr = BedStepCMB(user.m,&user); CHKERRQ(ierr);
          ierr = BedStepExactThickness(user.Hexact,&user); CHKERRQ(ierr);
      } else {
          SETERRQ(PETSC_COMM_WORLD,1,"ERROR: one of user.[dome,bedstep] must be TRUE since user.read is FALSE...\n");
      }
      ierr = ChopScaleSMBforInitial(user.Hinitial,&user); CHKERRQ(ierr);
  }

  if (user.showdata) {
      ierr = ShowFields(&user); CHKERRQ(ierr);
  }

  // setup local copy of bed
  ierr = DMCreateLocalVector(user.da,&user.bloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(user.da,user.b,INSERT_VALUES,user.bloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(user.da,user.b,INSERT_VALUES,user.bloc);CHKERRQ(ierr);

  // initialize the SNESVI
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetDM(snes,user.da);CHKERRQ(ierr);
  ierr = DMDASNESSetFunctionLocal(user.da,INSERT_VALUES,
                (DMDASNESFunction)FormFunctionLocal,&user); CHKERRQ(ierr);
  ierr = DMDASNESSetJacobianLocal(user.da,
                (DMDASNESJacobian)FormJacobianLocal,&user); CHKERRQ(ierr);
  ierr = SNESSetType(snes,SNESVINEWTONRSLS);CHKERRQ(ierr);
  ierr = SNESVISetComputeVariableBounds(snes,&FormBounds);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  if (user.mtrue)
      myPrintf(&user,    "solving by true Mahaffy method ...\n");
  else
      if (user.lambda == 0.0)
          myPrintf(&user,"solving by M* method _without_ upwinding ...\n");
      else
          myPrintf(&user,"solving by M* method ...\n");

  SNESConvergedReason reason;
  gettimeofday(&user.starttime, NULL);

  // time-stepping loop
  // note "Hinitial" is really "H^{l-1}", i.e. the value of the solution at the
  // last time step if time stepping
  for (l = 0; l < user.numBEsteps; l++) {
      ierr = VecCopy(user.Hinitial,H); CHKERRQ(ierr);
      // continuation loop
      for (m = user.cs.start; m<user.cs.end; m++) {
          user.eps = epsCS(m,&(user.cs));
          ierr = VecCopy(H,Htry); CHKERRQ(ierr);
          ierr = SNESAttempt(&snes,Htry,m,&reason,&user);CHKERRQ(ierr);
          if (reason < 0) {
              if (user.divergetryagain) {
                  myPrintf(&user,
                      "         turning on recovery mode (backward Euler time step of %.2f a) and trying again ...\n",
                      user.dtBE/user.secpera);
                  user.doBEsteps = PETSC_TRUE;
                  user.cs.goodm = m-1;
                  ierr = VecCopy(H,user.Hinitial); CHKERRQ(ierr);      
                  ierr = VecCopy(H,Htry); CHKERRQ(ierr);
                  ierr = SNESAttempt(&snes,Htry,m,&reason,&user);CHKERRQ(ierr);
              }
              if (reason < 0) {
                  if (m>0) // record last successful eps
                      user.eps = epsCS(m-1,&(user.cs));
                  break;
              }
          }
          if ((user.divergetryagain) && (user.doBEsteps))
              user.recoverycount++;
          else
              user.cs.goodm = m;
          ierr = VecCopy(Htry,H); CHKERRQ(ierr);
          ierr = StdoutReport(H,&user); CHKERRQ(ierr);
      }  // for m
      if (reason >= 0) {
          if (user.numBEsteps > 1)
              myPrintf(&user,"t = %.2f a: completed time step %d of %d\n",
                       (l+1)*user.dtBE/user.secpera,l+1,user.numBEsteps);
          ierr = VecCopy(H,user.Hinitial); CHKERRQ(ierr);
      } else
          break;
  } // for l

  gettimeofday(&user.endtime, NULL);

  if (user.history) {
      ierr = WriteHistoryFile(H,"history.txt",argc,argv,&user); CHKERRQ(ierr);
  }

  if (user.dump) {
      Vec r;
      ierr = SNESGetFunction(snes,&r,NULL,NULL); CHKERRQ(ierr);
      ierr = DumpToFile(H,r,&user); CHKERRQ(ierr);
  }

  if ((user.averr) || (user.maxerr)) {
      if (reason < 0)
          PetscPrintf(PETSC_COMM_WORLD,"%s\n",SNESConvergedReasons[reason]);
      else {
          PetscReal enorminf, enorm1;
          ierr = GetErrors(H, &user, &enorminf, &enorm1); CHKERRQ(ierr);
          if (user.averr)
              PetscPrintf(PETSC_COMM_WORLD,"%.14e\n",(double)enorm1 / (info.mx * info.my));
          if (user.maxerr)
              PetscPrintf(PETSC_COMM_WORLD,"%.14e\n",(double)enorminf);
      }
  }

  ierr = VecDestroy(&user.bloc);CHKERRQ(ierr);
  ierr = VecDestroy(&user.ds.Dnodemax);CHKERRQ(ierr);
  ierr = VecDestroy(&user.ds.Wmagnodemax);CHKERRQ(ierr);
  ierr = VecDestroy(&user.m);CHKERRQ(ierr);
  ierr = VecDestroy(&user.b);CHKERRQ(ierr);
  ierr = VecDestroy(&user.Hexact);CHKERRQ(ierr);
  ierr = VecDestroy(&user.Hinitial);CHKERRQ(ierr);
  ierr = VecDestroy(&Htry);CHKERRQ(ierr);
  ierr = VecDestroy(&H);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = DMDestroy(&user.quadda);CHKERRQ(ierr);
  ierr = DMDestroy(&user.sixteenda);CHKERRQ(ierr);
  ierr = DMDestroy(&user.da);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}


// set initial H by chop & scale SMB
PetscErrorCode ChopScaleSMBforInitial(Vec Hinitial, AppCtx *user) {
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = VecCopy(user->m,Hinitial); CHKERRQ(ierr);
  ierr = VecTrueChop(Hinitial,0.0); CHKERRQ(ierr);
  ierr = VecScale(Hinitial,user->initmagic * user->secpera); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


//  for call-back: tell SNESVI (variational inequality) that we want
//    0.0 <= H < +infinity
PetscErrorCode FormBounds(SNES snes, Vec Xl, Vec Xu) {
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = VecSet(Xl,0.0); CHKERRQ(ierr);
  ierr = VecSet(Xu,PETSC_INFINITY); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode ProcessOptions(AppCtx *user) {
  PetscErrorCode ierr;
  PetscBool      domechosen, dtBEset, notryset;
  char           histprefix[512];

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"mah_","options to mahaffy","");CHKERRQ(ierr);
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
      "-dtBE", "duration in years of time step used in Backward Euler recovery",
      NULL,user->dtBE/user->secpera,&user->dtBE,&dtBEset);CHKERRQ(ierr);
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
  ierr = PetscOptionsBool(
      "-readinitial", "use read initial H instead of generating guess in usual way",
      NULL,user->readinitial,&user->readinitial,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-showdata", "use PETSc X viewers to show b, m, and exact/observed H",
      NULL,user->showdata,&user->showdata,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-silent", "run silent (print nothing)",
      NULL,user->silent,&user->silent,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt(
      "-steps", "number of Backward Euler time-steps to take; default=1 is for steady-state",
      NULL,user->numBEsteps,&user->numBEsteps,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-swapxy", "swap coordinates x and y when building bedrock step exact solution",
      NULL,user->swapxy,&user->swapxy,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-true", "use true Mahaffy method, not default M*",
      NULL,user->mtrue,&user->mtrue,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  // enforce consistency of cases
  if (dtBEset)
      user->dtBE *= user->secpera;
  if (user->numBEsteps > 1)
      user->doBEsteps = PETSC_TRUE;
  if ((user->averr) || (user->maxerr))
      user->silent = PETSC_TRUE;
  if (user->read) {
      if (domechosen) {
        SETERRQ(PETSC_COMM_WORLD,1,"ERROR option conflict: both -mah_dome and -mah_read not allowed\n");
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
      user->lambda = 0.0;
  if (user->dump)
      user->history = PETSC_TRUE;
  else if (user->history)
      strcpy(user->figsprefix,histprefix);

  ierr = OptionsCS(&(user->cs)); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}




