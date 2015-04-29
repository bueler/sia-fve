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
   mpiexec -n 6 ./mahaffy -da_refine 4 -mah_dump teststeps/ -mah_notry
   # above diverges at 12 with          errors:  max =  268.69 m,  av =    8.04 m,  voldiff% =  0.70
   mpiexec -n 6 ./mahaffy -mah_read teststeps/unnamed.dat -mah_readinitial teststeps/unnamed.dat -cs_start 10 -mah_steps 100 -mah_dt 100.0 -mah_notry
   # above completes 10000 a run with   errors:  max =  217.21 m,  av =    0.74 m,  voldiff% =  0.00
*/

#include <math.h>
#include <petscdmda.h>
#include <petscsnes.h>

#include "base/appctx.h"
#include "base/continuationscheme.h"
#include "base/exactsia.h"
#include "base/io.h"
#include "base/solver.h"

extern PetscErrorCode ChopScaleSMBforInitial(Vec,AppCtx*);
extern PetscErrorCode FormBounds(SNES,Vec,Vec);
extern PetscErrorCode Step(Vec,SNES*,ContinuationScheme*,SNESConvergedReason*,AppCtx*);


int main(int argc,char **argv) {
  PetscErrorCode      ierr;
  SNES                snes;
  Vec                 H;
  AppCtx              user;
  ContinuationScheme  cs;
  DMDALocalInfo       info;
  PetscReal           dx,dy;
  SNESConvergedReason reason;

  PetscInitialize(&argc,&argv,(char*)0,help);

  ierr = SetFromOptionsAppCtx("mah_",&user); CHKERRQ(ierr);
  ierr = SetFromOptionsCS("cs_",&cs); CHKERRQ(ierr);
  user.cs = &cs;

  if (user.read) {
      myPrintf(&user,
          "reading dimensions and grid spacing from %s ...\n",user.readname);
      ierr = ReadDimensions(&user); CHKERRQ(ierr);  // sets user.[Nx,Ny,Lx,Ly,dx]
  } else {
      if (user.dome)
          DomeDefaultGrid(&user);
      else
          BedStepDefaultGrid(&user);
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
      myPrintf(&user,"reading b, m, Hexact (or Hobserved) from %s ...\n", user.readname);
      ierr = ReadDataVecs(&user); CHKERRQ(ierr);
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
  }

  // fill user.Hinitial according to either -mah_readinitial foo.dat or chop-scale-SMB
  if (user.readinitial) {
      myPrintf(&user,"  reading Hinitial from %s ...\n", user.readinitialname);
      ierr = ReadInitialHVec(&user); CHKERRQ(ierr);
  } else {
      myPrintf(&user,"  generating Hinitial by chop-and-scale of SMB ...\n");
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

  myPrintf(&user,    "solving %s by %s method ...\n",
           (user.dtres > 0.0) ? "backward-Euler time-step problems" : "steady-state problem",
           (user.mtrue) ? "true Mahaffy" :
               ((user.lambda == 0.0) ? "M*-without-upwinding" : "M*"));
  if (user.dtres > 0.0) {
      myPrintf(&user,"  time-stepping on interval [0.0 a,%.2f a] with initial step %.2f a\n",
               user.T/user.secpera,user.dtres/user.secpera);
  }

  ierr = VecCopy(user.Hinitial,H); CHKERRQ(ierr);

  gettimeofday(&user.starttime, NULL);

  if (user.dtres > 0.0) {  // time-stepping
      // note "Hinitial" is really "H^{l-1}", the solution at the last time step
      PetscReal  tcurrent = 0.0,
                 dtgoal   = user.dtres;
      while (tcurrent < user.T) {
          PetscInt reducecount = 10;
          user.dtres = PetscMin(user.T - tcurrent, dtgoal);
          user.dtjac = user.dtres;
          ierr = Step(H,&snes,&cs,&reason,&user); CHKERRQ(ierr);
          while (reason < 0) {
              if (reducecount == 0) {
                  SETERRQ(PETSC_COMM_WORLD,1,"ERROR:  time-step failure ... giving up on reducing time step\n");
              }
              reducecount--;
              ierr = VecCopy(user.Hinitial,H); CHKERRQ(ierr);
              user.dtres /= 2.0;
              user.dtjac = user.dtres;
              ierr = Step(H,&snes,&cs,&reason,&user); CHKERRQ(ierr);
          }
          ierr = VecCopy(H,user.Hinitial); CHKERRQ(ierr);
          tcurrent += user.dtres;
          myPrintf(&user,"t = %.2f a: completed time step of duration %.2f a in interval [0.0 a,%.2f a]\n",
                   tcurrent/user.secpera,user.dtres/user.secpera,user.T/user.secpera);
      }
  } else { // steady-state
      ierr = Step(H,&snes,&cs,&reason,&user); CHKERRQ(ierr);
      if (reason >= 0) {
          ierr = VecCopy(H,user.Hinitial); CHKERRQ(ierr);
      }
  }

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

// solve the steady-state problem or do one time step:
//   * applies the continuation scheme
//   * applies recovery if in steady-state and if desired
PetscErrorCode Step(Vec H, SNES *snes, ContinuationScheme *cs, SNESConvergedReason *reason, AppCtx *user) {
  PetscErrorCode ierr;
  Vec            Htry;
  PetscInt       m;

  PetscFunctionBeginUser;
  ierr = VecDuplicate(H,&Htry);CHKERRQ(ierr);

  for (m = startCS(cs); m < endCS(cs); m++) {
      user->eps = epsCS(m,cs);
      ierr = VecCopy(H,Htry); CHKERRQ(ierr);
      ierr = SNESAttempt(snes,Htry,m,reason,user);CHKERRQ(ierr);
      if (*reason < 0) {
          if ((user->divergetryagain) && (user->recoverycount == 0) && (user->dtres <= 0.0)) {
              myPrintf(user,
                  "         turning on steady-state recovery mode (backward Euler time step of %.2f a) ...\n",
                  user->dtrecovery/user->secpera);
              user->dtres = user->dtrecovery;
              user->dtjac = user->dtrecovery;
              user->recoverycount = 1;
              user->goodm = m-1;
              myPrintf(user,
                  "         trying again ...\n");
              ierr = VecCopy(H,user->Hinitial); CHKERRQ(ierr);
              ierr = VecCopy(H,Htry); CHKERRQ(ierr);
              ierr = SNESAttempt(snes,Htry,m,reason,user);CHKERRQ(ierr);
          }
          if (*reason < 0) {
              if (m>0) // record last successful eps
                  user->eps = epsCS(m-1,cs);
              break;
          }
      } else if (user->recoverycount > 0)
          user->recoverycount++;
      // actions when successful:
      ierr = VecCopy(Htry,H); CHKERRQ(ierr);
      ierr = StdoutReport(H,user); CHKERRQ(ierr);
      if (user->recoverycount == 0)
          user->goodm = m;
  }

  ierr = VecDestroy(&Htry);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
