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
   for LEV in 0 1 2 3 4; do  mpiexec -n 6 ./mahaffy -da_refine $LEV -snes_type vinewtonssls -snes_max_it 200 -pc_type asm -sub_pc_type lu; done

Successes:
   mpiexec -n 6 ./mahaffy -pc_type mg -da_refine 5 -snes_max_it 200 -snes_monitor  # DIVERGED_LINE_SEARCH at 12
   mpiexec -n 6 ./mahaffy -da_refine 6 -pc_type asm -sub_pc_type lu -snes_max_it 200

Divergence:
   ./mahaffy -da_refine 0                   # no divergence
   ./mahaffy -da_refine 1                   # no divergence
   ./mahaffy -da_refine 2                   # divergence at 12
   ./mahaffy -da_refine 3                   # divergence at 11
   ./mahaffy -da_refine 4 -snes_max_it 200  # divergence at 11
   ./mahaffy -mah_bedstep -mah_D0 0.01 -da_refine 0    # DIVERGED_LINE_SEARCH at 12
   ./mahaffy -mah_bedstep -mah_D0 0.01 -da_refine 1 -snes_max_it 1000  # DIVERGED_FUNCTION_COUNT at 13
   ./mahaffy -mah_bedstep -mah_D0 0.01 -da_refine 2    # DIVERGED_LINE_SEARCH at 13
   ./mahaffy -mah_bedstep -mah_D0 0.01 -da_refine 3    # DIVERGED_LINEAR_SOLVE at 4
   ./mahaffy -mah_bedstep -mah_D0 0.01 -da_refine 4    # DIVERGED_LINEAR_SOLVE at 4
   # zero pivot at 13:
   mpiexec -n 6 ./mahaffy -mah_bedstep -mah_D0 0.01 -da_refine 4 -snes_max_it 1000 -pc_type asm -sub_pc_type lu
*/

#include <math.h>
#include <petscdmda.h>
#include <petscsnes.h>

#include "base/mahaffyctx.h"
#include "base/continuationscheme.h"
#include "base/exactsia.h"
#include "base/io.h"
#include "base/q1op.h"
#include "base/sia.h"

extern PetscErrorCode FormBounds(SNES,Vec,Vec);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*,PetscScalar**,PetscScalar**,AppCtx*);
extern PetscErrorCode FormJacobianLocal(DMDALocalInfo*,PetscScalar**,Mat,Mat,AppCtx*);
extern PetscErrorCode SNESAttempt(SNES*,Vec,PetscBool,PetscInt,SNESConvergedReason*,AppCtx*);
extern PetscErrorCode ProcessOptions(AppCtx*);
extern PetscErrorCode StateReport(Vec,DMDALocalInfo*,AppCtx*);
extern PetscErrorCode ChopScaleSMBforInitial(Vec,AppCtx*);
extern void fillscheds(AppCtx*);

// indexing of the 8 quadrature points along the boundary of the control volume in M*
// point s=0,...,7 is in element (j,k) = (j+je[s],k+ke[s])
static const PetscInt  je[8] = {0,  0, -1, -1, -1, -1,  0,  0},
                       ke[8] = {0,  0,  0,  0, -1, -1, -1, -1},
                       ce[8] = {0,  3,  1,  0,  2,  1,  3,  2};

// coefficients of quadrature evaluations along the boundary of the control volume in M*
#define FLUXINTCOEFFS \
const PetscReal coeff[8] = {dy/2, dx/2, dx/2, -dy/2, -dy/2, -dx/2, -dx/2, dy/2};

// direction of flux at 4 points in each element
static const PetscBool xdire[4] = {PETSC_TRUE, PETSC_FALSE, PETSC_TRUE, PETSC_FALSE};

// local (element-wise) coords of quadrature points for M*
static const PetscReal locx[4] = {  0.5, 0.75,  0.5, 0.25},
                       locy[4] = { 0.25,  0.5, 0.75,  0.5};

// quadrature points for classical (true) Mahaffy; these are not used in Jacobian
static const PetscReal locxtrue[4] = { 0.5, 1.0, 0.5, 0.0},
                       locytrue[4] = { 0.0, 0.5, 1.0, 0.5},
                       locxnbr[4]  = { 0.5, 0.0, 0.5, 1.0},
                       locynbr[4]  = { 1.0, 0.5, 0.0, 0.5};
static const PetscInt  jnbr[4] = { 0,  1,  0, -1},
                       knbr[4] = {-1,  0,  1,  0};


int main(int argc,char **argv) {
  PetscErrorCode      ierr;
  SNES                snes;
  Vec                 H, Htry;
  AppCtx              user;
  DMDALocalInfo       info;
  PetscReal           dx,dy;
  PetscInt            m;

  PetscInitialize(&argc,&argv,(char*)0,help);

  // default settings of parameters
  user.n      = 3.0;
  user.g      = 9.81;       // m/s^2
  user.rho    = 910.0;      // kg/m^3
  user.secpera= 31556926.0;
  user.A      = 1.0e-16/user.secpera; // = 3.17e-24  1/(Pa^3 s); EISMINT I value

  user.initmagic = 1000.0;  // a
  user.delta  = 1.0e-4;

  user.lambda = 0.25;  // amount of upwinding; some trial-and-error with bedstep soln; 0.1 gives some Newton convergence problem on refined grid (=125m) but this does not; earlier M* was 0.5 here
  user.dtBE   = 1.0 * user.secpera;  // default 1 year time step for Backward Euler; used only in recovery

  user.mtrue      = PETSC_FALSE;

  user.dome       = PETSC_TRUE;  // defaults to this case
  user.bedstep    = PETSC_FALSE;
  user.swapxy     = PETSC_FALSE;
  user.divergetryagain = PETSC_TRUE;
  user.dorecovery      = PETSC_FALSE;
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
  ierr = VecCopy(user.Hinitial,H); CHKERRQ(ierr);
  for (m = user.cs.start; m<user.cs.end; m++) {
      user.eps = epsCS(m,&(user.cs));
      ierr = VecCopy(H,Htry); CHKERRQ(ierr);
      ierr = SNESAttempt(&snes,Htry,PETSC_FALSE,m,&reason,&user);CHKERRQ(ierr);
      if (reason < 0) {
          if ((user.divergetryagain) && (!user.dorecovery)) {
              myPrintf(&user,
                  "         turning on recovery mode (backward Euler step of %.2f a) and trying again ...\n",
                  user.dtBE/user.secpera);
              user.dorecovery = PETSC_TRUE;
              ierr = VecCopy(H,user.Hinitial); CHKERRQ(ierr);  // only in case we need to recover
              ierr = VecCopy(H,Htry); CHKERRQ(ierr);
              ierr = SNESAttempt(&snes,Htry,PETSC_TRUE,m,&reason,&user);CHKERRQ(ierr);
          }
          if (reason < 0) {
              if (m>0) // record last successful eps
                  user.eps = epsCS(m-1,&(user.cs));
              break;
          }
      }
      ierr = VecCopy(Htry,H); CHKERRQ(ierr);
      ierr = StdoutReport(H,&user); CHKERRQ(ierr);
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


/* Loop over locally-owned elements, including ghosts, checking
   nonnegativity of thickness.  Stops with error if not.  */
PetscErrorCode checkadmissible(DMDALocalInfo *info, PetscScalar **H) {
  PetscInt        j, k;
  PetscFunctionBeginUser;
  for (k = info->ys-1; k < info->ys + info->ym + 1; k++) {
      for (j = info->xs-1; j < info->xs + info->xm + 1; j++) {
          if (H[k][j] < 0.0)
              SETERRQ3(PETSC_COMM_WORLD,1,"ERROR: inadmissible H[%d][%d] = %.3e < 0\n",k,j,H[k][j]);
      }
  }
  PetscFunctionReturn(0);
}


// averages two gradients; only used in true Mahaffy
Grad gradav(Grad g1, Grad g2) {
  Grad gav;
  gav.x = (g1.x + g2.x) / 2.0;
  gav.y = (g1.y + g2.y) / 2.0;
  return gav;
}


/* For call-back by SNES using DMDA info.

Evaluates residual FF on local process patch:
   FF_{j,k} = \int_{\partial V_{j,k}} \mathbf{q} \cdot \mathbf{n} - m_{j,k} \Delta x \Delta y
where V_{j,k} is the control volume centered at (x_j,y_k).

Regarding indexing the location along the boundary of the control volume where
flux is evaluated, this shows four elements and one control volume centered
at (x_j,y_k).  The boundary of the control volume has 8 points, numbered s=0,...,7:
   -------------------
  |         |         |
  |    ..2..|..1..    |
  |   3:    |    :0   |
k |--------- ---------|
  |   4:    |    :7   |
  |    ..5..|..6..    |
  |         |         |
   -------------------
            j

Regarding flux-component indexing on the element indexed by (j,k) node, as shown,
the value  aq[k][j][c]  for c=0,1,2,3, is an x-component at "*" and a y-component
at "%":
   -------------------
  |         :         |
  |         *2        |
  |    3    :    1    |
  |....%.... ....%....|
  |         :         |
  |         *0        |
  |         :         |
  @-------------------
(j,k)
*/
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info, PetscScalar **aH, PetscScalar **FF, AppCtx *user) {
  PetscErrorCode  ierr;
  const PetscReal dx = user->dx, dy = dx;
  FLUXINTCOEFFS
  const PetscBool upwind = (user->lambda > 0.0);
  const PetscReal upmin = (1.0 - user->lambda) * 0.5,
                  upmax = (1.0 + user->lambda) * 0.5;
  PetscInt        j, k;
  PetscReal       **am, **ab, **aHprev, **aDnodemax, **aWmagnodemax,
                  ***aqquad,  ***aDquad, ***aWquad;
  Vec             qloc, Dloc, Wloc;
  DiagnosticScheme *ds = &(user->ds);

  PetscFunctionBeginUser;
  if (user->checkadmissible) {
      ierr = checkadmissible(info,aH); CHKERRQ(ierr);
  }

  if (!user->nodiag) {
      ds->avD      = 0.0;
      ds->avDcount = 0;
      ds->maxD     = 0.0;
      ierr = DMGetLocalVector(user->quadda,&Dloc);CHKERRQ(ierr);
      ierr = DMGetLocalVector(user->quadda,&Wloc);CHKERRQ(ierr);
      ierr = DMDAVecGetArray(user->da, ds->Dnodemax, &aDnodemax);CHKERRQ(ierr);
      ierr = DMDAVecGetArray(user->da, ds->Wmagnodemax, &aWmagnodemax);CHKERRQ(ierr);
      ierr = DMDAVecGetArrayDOF(user->quadda, Dloc, &aDquad);CHKERRQ(ierr);
      ierr = DMDAVecGetArrayDOF(user->quadda, Wloc, &aWquad);CHKERRQ(ierr);
  }

  // need stencil width on locally-computed q
  ierr = DMGetLocalVector(user->quadda,&qloc);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->bloc, &ab);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->m, &am);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->Hinitial, &aHprev);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(user->quadda, qloc, &aqquad);CHKERRQ(ierr);
  // loop over locally-owned elements, including ghosts, to get fluxes at
  // c = 0,1,2,3 points in element;  note start at (xs-1,ys-1)
  for (k = info->ys-1; k < info->ys + info->ym; k++) {
      for (j = info->xs-1; j < info->xs + info->xm; j++) {
          PetscInt  c;
          PetscReal H, Hup;
          Grad      gH, gb;
          for (c=0; c<4; c++) {
              if (user->mtrue) {  // true Mahaffy method
                  // this implementation is at least a factor of two inefficient
                  H = fieldatpt(j,k,locxtrue[c],locytrue[c],aH);
                  gH = gradfatpt(j,k,locxtrue[c],locytrue[c],dx,dy,aH);
                  gb = gradfatpt(j,k,locxtrue[c],locytrue[c],dx,dy,ab);
                  gH = gradav(gH,gradfatpt(j+jnbr[c],k+knbr[c],locxnbr[c],locynbr[c],dx,dy,aH));
                  gb = gradav(gb,gradfatpt(j+jnbr[c],k+knbr[c],locxnbr[c],locynbr[c],dx,dy,ab));
                  Hup = H; // no upwinding allowed in true Mahaffy
              } else { // default M* method
                  H  = fieldatpt(j,k,locx[c],locy[c],aH);
                  gH = gradfatpt(j,k,locx[c],locy[c],dx,dy,aH);
                  gb = gradfatpt(j,k,locx[c],locy[c],dx,dy,ab);
                  if (upwind) {
                      PetscReal lxup = locx[c], lyup = locy[c];
                      if (xdire[c] == PETSC_TRUE)
                          lxup = (gb.x <= 0.0) ? upmin : upmax;
                      else
                          lyup = (gb.y <= 0.0) ? upmin : upmax;
                      Hup = fieldatpt(j,k,lxup,lyup,aH);
                  } else
                      Hup = H;
              }
              if (user->nodiag)
                  aqquad[k][j][c] = getflux(gH,gb,H,Hup,xdire[c],user);
              else
                  aqquad[k][j][c] = getfluxDIAGNOSTIC(gH,gb,H,Hup,xdire[c],user,
                                                      &(aDquad[k][j][c]),&(aWquad[k][j][c]));
          }
      }
  }
  // loop over nodes, not including ghosts, to get residual from quadature over
  // s = 0,1,...,7 points on boundary of control volume (rectangle) around node
  for (k=info->ys; k<info->ys+info->ym; k++) {
      for (j=info->xs; j<info->xs+info->xm; j++) {
          PetscInt s;
          // This is the integral over the control volume boundary using two
          // quadrature points on each side of of the four sides of the
          // rectangular control volume.
          // For M*: two instances of midpoint rule on each side, with two
          //         different values of aq[][] per side
          // For true Mahaffy: ditto, but the two values of aq[][] per side are
          //                   actually the same.
          FF[k][j] = - am[k][j] * dx * dy;
          for (s=0; s<8; s++)
              FF[k][j] += coeff[s] * aqquad[k+ke[s]][j+je[s]][ce[s]];
          if (user->dorecovery)
              FF[k][j] += (aH[k][j] - aHprev[k][j]) * dx * dy / user->dtBE;
          if (!user->nodiag) {
              // update diagnostics associated to diffusivity
              PetscReal Dmax = 0.0, Wmagmax = 0.0;
              for (s=0; s<8; s++) {
                  const PetscReal D = aDquad[k+ke[s]][j+je[s]][ce[s]];
                  ds->maxD        = PetscMax(ds->maxD, D);
                  Dmax            = PetscMax(Dmax, D);
                  ds->avD        += D;
                  ds->avDcount   += 1;
                  Wmagmax         = PetscMax(Wmagmax, aWquad[k+ke[s]][j+je[s]][ce[s]]);
              }
              aDnodemax[k][j] = Dmax;
              aWmagnodemax[k][j] = Wmagmax;
          }
      }
  }
  ierr = DMDAVecRestoreArray(user->da, user->m, &am);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, user->Hinitial, &aHprev);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, user->bloc, &ab);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(user->quadda, qloc, &aqquad);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(user->quadda,&qloc);CHKERRQ(ierr);

  if (!user->nodiag) {
      ierr = DMDAVecRestoreArray(user->da, ds->Dnodemax, &aDnodemax);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArray(user->da, ds->Wmagnodemax, &aWmagnodemax);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArrayDOF(user->quadda, Dloc, &aDquad);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArrayDOF(user->quadda, Wloc, &aWquad);CHKERRQ(ierr);
      ierr = DMRestoreLocalVector(user->quadda,&Dloc);CHKERRQ(ierr);
      ierr = DMRestoreLocalVector(user->quadda,&Wloc);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}


// allows consistent use of "j" and "k" for x,y directions, both in loops and MatSetValuesStencil
typedef struct {
  PetscInt foo,k,j,bar;
} MyStencil;

/* For call-back by SNES using DMDA info.

Evaluates Jacobian matrix on local process patch.

For examples see $PETSC_DIR/src/snes/examples/tutorials/ex5.c or ex9.c.
*/
PetscErrorCode FormJacobianLocal(DMDALocalInfo *info, PetscScalar **aH, Mat jac, Mat jacpre, AppCtx *user) {
  PetscErrorCode  ierr;
  const PetscReal dx = user->dx, dy = dx;
  FLUXINTCOEFFS
  const PetscBool upwind = (user->lambda > 0.0);
  const PetscReal upmin = (1.0 - user->lambda) * 0.5,
                  upmax = (1.0 + user->lambda) * 0.5;
  PetscInt        j, k;
  PetscReal       **ab, ***adQ;
  Vec             dQloc;
  MyStencil       col[33],row;
  PetscReal       val[33];

  PetscFunctionBeginUser;
  if (user->mtrue) {
      SETERRQ(PETSC_COMM_WORLD,1,"ERROR: analytical jacobian not ready in this cases ...\n"); }
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"[inside FormJacobianLocal()]\n"); CHKERRQ(ierr);
  if (user->checkadmissible) {
      ierr = checkadmissible(info,aH); CHKERRQ(ierr); }

  ierr = MatZeroEntries(jac); CHKERRQ(ierr);  // because using ADD_VALUES below

  ierr = DMGetLocalVector(user->sixteenda,&dQloc);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(user->da, user->bloc, &ab);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(user->sixteenda, dQloc, &adQ);CHKERRQ(ierr);
  // loop over locally-owned elements, including ghosts, to get DfluxDl for
  // l=0,1,2,3 at c = 0,1,2,3 points in element;  note start at (xs-1,ys-1)
  for (k = info->ys-1; k < info->ys + info->ym; k++) {
      for (j = info->xs-1; j < info->xs + info->xm; j++) {
          PetscInt  c, l;
          Grad      gH, gb;
          PetscReal H, Hup;
          for (c=0; c<4; c++) {
              PetscReal lxup = locx[c], lyup = locy[c];
              H  = fieldatpt(j,k,locx[c],locy[c],aH);
              gH = gradfatpt(j,k,locx[c],locy[c],dx,dy,aH);
              gb = gradfatpt(j,k,locx[c],locy[c],dx,dy,ab);
              Hup = H;
              if (upwind) {
                  if (xdire[c])
                      lxup = (gb.x <= 0.0) ? upmin : upmax;
                  else
                      lyup = (gb.y <= 0.0) ? upmin : upmax;
                  Hup = fieldatpt(j,k,lxup,lyup,aH);
              }
              for (l=0; l<4; l++) {
                  Grad      dgHdl;
                  PetscReal dHdl, dHupdl;
                  dgHdl  = dgradfatpt(l,j,k,locx[c],locy[c],dx,dy);
                  dHdl   = dfieldatpt(l,j,k,locx[c],locy[c]);
                  dHupdl = (upwind) ? dfieldatpt(l,j,k,lxup,lyup) : dHdl;
                  adQ[k][j][4*c+l] = DfluxDl(gH,gb,dgHdl,H,dHdl,Hup,dHupdl,xdire[c],user);
              }
          }
      }
  }
  // loop over nodes, not including ghosts, to get derivative of residual with respect to nodal value
  for (k=info->ys; k<info->ys+info->ym; k++) {
      row.k = k;
      for (j=info->xs; j<info->xs+info->xm; j++) {
          row.j = j;
          PetscInt s, u, v, l;
          for (s=0; s<8; s++) {
              u = j + je[s];
              v = k + ke[s];
              for (l=0; l<4; l++) {
                  const PetscInt djfroml[4] = { 0,  1,  1,  0},
                                 dkfroml[4] = { 0,  0,  1,  1};
                  col[4*s+l].j = u + djfroml[l];
                  col[4*s+l].k = v + dkfroml[l];
                  val[4*s+l] = coeff[s] * adQ[v][u][4*ce[s]+l];
              }
          }
          if (user->dorecovery) {
              // add another stencil for diagonal
              col[32].j = j;
              col[32].k = k;
              val[32]   = dx * dy / user->dtBE;
              ierr = MatSetValuesStencil(jac,1,(MatStencil*)&row,33,(MatStencil*)col,val,ADD_VALUES);CHKERRQ(ierr);
          } else {
              ierr = MatSetValuesStencil(jac,1,(MatStencil*)&row,32,(MatStencil*)col,val,ADD_VALUES);CHKERRQ(ierr);
          }
      }
  }
  ierr = DMDAVecRestoreArray(user->da, user->bloc, &ab);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(user->sixteenda, dQloc, &adQ);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(user->sixteenda,&dQloc);CHKERRQ(ierr);

  // Assemble matrix, using the 2-step process:
  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (jacpre != jac) {
    ierr = MatAssemblyBegin(jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }

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
  ierr = PetscOptionsBool(
      "-swapxy", "swap coordinates x and y when building bedrock step exact solution",
      NULL,user->swapxy,&user->swapxy,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-true", "use true Mahaffy method, not default M*",
      NULL,user->mtrue,&user->mtrue,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  // enforce consistency of cases
  if (dtBEset) {
      if (!user->divergetryagain) {
        SETERRQ(PETSC_COMM_WORLD,1,"ERROR setting -mah_dtBE has no effect if -mah_notry is also set\n");
      } else
        user->dtBE *= user->secpera;
  }
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

PetscErrorCode PetscIgnoreZEROPIVOTErrorHandler(MPI_Comm comm,int line,const char *fun,
                   const char *file,PetscErrorCode n,PetscErrorType p,const char *mess,void *ctx) {
   if ((n == PETSC_ERR_MAT_LU_ZRPVT) || (p == PETSC_ERR_MAT_LU_ZRPVT)) {
      int rank;
      MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
      //PetscPrintf(PETSC_COMM_SELF,"---- WARNING: intercepted PETSC_ERR_MAT_LU_ZRPVT on rank %d ----\n",rank);
      AppCtx* user = (AppCtx*)ctx;
      user->luzeropvterr = 1;
      return 0;
   } else {
      return PetscTraceBackErrorHandler(comm,line,fun,file,n,p,mess,ctx);
   }
}

// try a SNES solve and report on the result; H is modified
PetscErrorCode SNESAttempt(SNES *s, Vec H, PetscBool again, PetscInt m,
                           SNESConvergedReason *reason, AppCtx *user) {
  PetscErrorCode ierr;
  KSP            ksp;
  PetscInt       its, kspits, luzeropvterr;
  const PetscInt lureason = -99;
  const char     lureasons[30] = "DIVERGED_LU_ZERO_PIVOT";
  char           reasonstr[30];

  PetscFunctionBeginUser;
  user->luzeropvterr = 0;
  PetscPushErrorHandler(PetscIgnoreZEROPIVOTErrorHandler,user);
  SNESSolve(*s, NULL, H);
  PetscPopErrorHandler();
  ierr = MPI_Allreduce(&user->luzeropvterr,&luzeropvterr,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD); CHKERRQ(ierr);
  if (luzeropvterr > 0) {
      *reason = lureason;
  } else {
      ierr = SNESGetConvergedReason(*s,reason);CHKERRQ(ierr);
  }
  ierr = SNESGetIterationNumber(*s,&its);CHKERRQ(ierr);
  ierr = SNESGetKSP(*s,&ksp); CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp,&kspits); CHKERRQ(ierr);
  strcpy(reasonstr,(*reason==lureason) ? lureasons : SNESConvergedReasons[*reason]);
  if (again)
      myPrintf(user,"     %s again w. ",reasonstr);
  else
      myPrintf(user,"%3d. %s   with   ",m,reasonstr);
  myPrintf(user,"eps=%.2e ... %3d KSP (last) iters and %3d Newton iters%s\n",
           user->eps,kspits,its,(user->dorecovery) ? " (recovery mode)" : "");
  PetscFunctionReturn(0);
}

