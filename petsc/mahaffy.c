/* (C) 2015 Ed Bueler */

static const char help[] =
"Solves steady ice sheet problem in 2d:\n"
"    div (q^x,q^y) = m,\n"
"    (q^x,q^y) = - Gamma H^{n+2} |grad s|^{n-1} grad s,\n"
"where  H(x,y)  is the ice thickness,  b(x,y)  is the bed elevation,\n"
"and  s(x,y) = H(x,y) + b(x,y)  is surface elevation.\n"
"Note  n > 1  and  Gamma = 2 A (rho g)^n / (n+2).\n"
"Domain is  -Lx < x < Lx,  -Ly < y < Ly,  with periodic boundary conditions.\n\n"
"Computed by Q1 FVE method with FD evaluation of Jacobian (i.e. no analytical yet).\n"
"Compares M* improved scheme (default) and true Mahaffy schemes.\n"
"Uses SNESVI (-snes_type vinewtonrsls) with constraint  H(x,y) >= 0.\n\n"
"Three problem cases:\n"
"  (1) flat bed case where analytical solution is known\n"
"  (2) bedrock-step case where analytical solution is known\n"
"  (3) -mah_read foo.dat reads b and m from PETSc binary file; see grn/README.md.\n\n";

/* Usage help:
   ./mahaffy -help |grep mah_

Use one of these problem cases:
   ./mahaffy -mah_dome         # default problem
   ./mahaffy -mah_bedstep
   ./mahaffy -mah_read foo.dat # see README.md for Greenland example

Solution options:
   ./mahaffy -mah_true -da_grid_x 15 -da_grid_y 15
                               # use true Mahaffy (sans quadrature & upwind improvements), but with 5|mx and 5|my
   ./mahaffy -mah_noupwind     # do not use upwinding on  grad b  part of flux
   ./mahaffy -mah_Neps 4       # don't go all the way on continuation
   ./mahaffy -mah_D0 10.0      # change constant diffusivity, used in continuation, to other value than
                               # default = 1.0;  big may be good for convergence, esp. w upwinding
   ./mahaffy -mah_divergetryagain # on SNES diverge, try again with eps *= 1.5

PETSc solver variations:
   ./mahaffy -mah_forceadmissible -snes_type vinewtonssls

Feedback on solution process:
   ./mahaffy -da_refine 1 -snes_vi_monitor  # widen screen to see SNESVI monitor output
   ./mahaffy -da_refine 3 -snes_monitor -snes_monitor_solution -snes_monitor_residual -draw_pause 1

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
   mpiexec -n 6 ./mahaffy -mah_bedstep -mah_D0 0.01 -da_refine 4 \  # zero pivot at 13
       -snes_max_it 1000 -pc_type asm -sub_pc_type lu

Generate .png figs:
   mkdir foo/
   ./mahaffy -mah_dump foo/
   cd foo/
   python ../figsmahaffy.py
*/

#include <math.h>
#include <petscdmda.h>
#include <petscsnes.h>

#include "mahaffyctx.h"
#include "exactsia.h"
#include "io.h"
#include "q1op.h"
#include "sia.h"

extern PetscErrorCode FormBounds(SNES,Vec,Vec);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*,PetscScalar**,PetscScalar**,AppCtx*);
// extern PetscErrorCode FormJacobianLocal(DMDALocalInfo*,PetscScalar**,Mat,Mat,AppCtx*); // see SNES ex9.c

extern PetscErrorCode SNESboot(SNES*,AppCtx*);
extern PetscErrorCode SNESAttempt(SNES,Vec,PetscInt*,SNESConvergedReason*);

extern PetscErrorCode ExplicitStepSmoother(Vec,AppCtx*);

extern PetscErrorCode ProcessOptions(AppCtx*);
extern PetscErrorCode StateReport(Vec,DMDALocalInfo*,AppCtx*);

// indexing of the 8 points along the boundary of the control volume
// point s=0,...,7 is in element (j,k) = (j+je[s],k+ke[s])
static const PetscInt  je[8] = {0,  0, -1, -1, -1, -1,  0,  0},
                       ke[8] = {0,  0,  0,  0, -1, -1, -1, -1},
                       ce[8] = {0,  3,  1,  0,  2,  1,  3,  2};
// direction of flux at 4 points in each element
static const PetscBool xdire[4] = {PETSC_TRUE, PETSC_FALSE, PETSC_TRUE, PETSC_FALSE};

int main(int argc,char **argv) {
  PetscErrorCode      ierr;
  SNES                snes;
  Vec                 H, Htry;
  AppCtx              user;
  DMDALocalInfo       info;
  PetscReal           dx,dy;

  PetscInitialize(&argc,&argv,(char*)0,help);

  // default settings of parameters
  user.n      = 3.0;
  user.g      = 9.81;       // m/s^2
  user.rho    = 910.0;      // kg/m^3
  user.secpera= 31556926.0;
  user.A      = 1.0e-16/user.secpera; // = 3.17e-24  1/(Pa^3 s); EISMINT I value
  user.Gamma  = 2.0 * PetscPowReal(user.rho*user.g,user.n) * user.A / (user.n+2.0);

  user.eps    = 0.0;
  user.D0     = 1.0;        // m^2 / s
  user.Neps   = 13;

  user.mtrue      = PETSC_FALSE;
  user.noupwind   = PETSC_FALSE;
  user.upwindfull = PETSC_FALSE;

  user.dome       = PETSC_TRUE;  // defaults to this case
  user.bedstep    = PETSC_FALSE;
  user.swapxy     = PETSC_FALSE;
  user.divergetryagain = PETSC_FALSE;
  user.checkadmissible = PETSC_FALSE;
  user.forceadmissible = PETSC_FALSE;

  user.read       = PETSC_FALSE;
  user.showdata   = PETSC_FALSE;
  user.history    = PETSC_FALSE;
  user.dump       = PETSC_FALSE;

  strcpy(user.figsprefix,"PREFIX/");  // dummies improve "mahaffy -help" appearance
  strcpy(user.readname,"FILENAME");

  ierr = ProcessOptions(&user); CHKERRQ(ierr);

  if (user.read == PETSC_TRUE) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
          "reading dimensions and grid spacing from %s ...\n",user.readname); CHKERRQ(ierr);
      ierr = ReadDimensions(&user); CHKERRQ(ierr);  // sets user.[Nx,Ny,Lx,Ly,dx]
  } else {
      if (user.dome == PETSC_TRUE) {
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
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      "solving on [-Lx,Lx]x[-Ly,Ly] with  Lx=%.3f km  and  Ly=%.3f km\n"
      "grid of  %d x %d  points with spacing  dx = %.6f km ...\n",
      user.Lx/1000.0,user.Ly/1000.0,info.mx,info.my,dx/1000.0); CHKERRQ(ierr);

  // this DMDA is used for evaluating flux components at quad points on elements
  ierr = DMDACreate2d(PETSC_COMM_WORLD,
                      DM_BOUNDARY_PERIODIC,DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX,
                      info.mx,info.my,PETSC_DECIDE,PETSC_DECIDE,
                      4, (user.mtrue==PETSC_TRUE) ? 2 : 1,  // dof=1,  stencilwidth=1
                      NULL,NULL,&user.quadda);
  ierr = DMSetApplicationContext(user.quadda, &user);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(user.da,&H);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)H,"thickness solution H"); CHKERRQ(ierr);

  ierr = VecDuplicate(H,&Htry); CHKERRQ(ierr);

  ierr = VecDuplicate(H,&user.b); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)(user.b),"bed elevation b"); CHKERRQ(ierr);
  ierr = VecDuplicate(H,&user.m); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)(user.m),"surface mass balance m"); CHKERRQ(ierr);
  ierr = VecDuplicate(H,&user.Hexact); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)H,"exact/observed thickness H"); CHKERRQ(ierr);

  // fill user.[b,m,Hexact] according to 3 choices: data, dome exact soln, JSA exact soln
  if (user.read == PETSC_TRUE) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
          "reading b, m, Hexact from %s ...\n", user.readname); CHKERRQ(ierr);
      ierr = ReadDataVecs(&user); CHKERRQ(ierr);
  } else {
      if (user.dome == PETSC_TRUE) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,
             "generating b, m, Hexact from dome formulas in Bueler (2003) ...\n"); CHKERRQ(ierr);
          ierr = VecSet(user.b,0.0); CHKERRQ(ierr);
          ierr = DomeCMB(user.m,&user); CHKERRQ(ierr);
          ierr = DomeExactThickness(user.Hexact,&user); CHKERRQ(ierr);
      } else if (user.bedstep == PETSC_TRUE) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,
             "generating b, m, Hexact from bedrock step formulas in Jarosch et al (2013) ...\n"); CHKERRQ(ierr);
          ierr = BedStepBed(user.b,&user); CHKERRQ(ierr);
          ierr = BedStepCMB(user.m,&user); CHKERRQ(ierr);
          ierr = BedStepExactThickness(user.Hexact,&user); CHKERRQ(ierr);
      } else {
          SETERRQ(PETSC_COMM_WORLD,1,"ERROR: one of user.[dome,bedstep] must be TRUE since user.read is FALSE...\n");
      }
  }

  if (user.showdata == PETSC_TRUE) {
      ierr = ShowFields(&user); CHKERRQ(ierr);
  }

  // initialize by chop & scale SMB
  ierr = VecCopy(user.m,H); CHKERRQ(ierr);
  ierr = VecChop(H,0.0); CHKERRQ(ierr);
  ierr = VecScale(H,1000.0*user.secpera); CHKERRQ(ierr);  // FIXME make user.initializemagic
  // alternatives:
  //ierr = [VERIF]ExactThickness(H,&user);CHKERRQ(ierr);
  //ierr = VecSet(H,0.0); CHKERRQ(ierr);

  if (user.mtrue == PETSC_TRUE) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,    "solving by true Mahaffy method ...\n"); CHKERRQ(ierr);
  } else {
      if (user.noupwind == PETSC_TRUE) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"solving by M* method _without_ upwinding ...\n"); CHKERRQ(ierr);
      } else {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"solving by M* method ...\n"); CHKERRQ(ierr);
      }
  }

  KSP        ksp;
  PetscInt   its, kspits, m;
  PetscReal  eps_sched[13] = {1.0,    0.5,    0.2,     0.1,     0.05,
                              0.02,   0.01,   0.005,   0.002,   0.001,
                              0.0005, 0.0002, 0.0000};
  SNESConvergedReason reason;
  ierr = SNESboot(&snes,&user); CHKERRQ(ierr);
  gettimeofday(&user.starttime, NULL);
  for (m = 0; m<user.Neps; m++) {
      user.eps = eps_sched[m];
      ierr = VecCopy(H,Htry); CHKERRQ(ierr);
      ierr = SNESAttempt(snes,Htry,&its,&reason);CHKERRQ(ierr);
      ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
      ierr = KSPGetIterationNumber(ksp,&kspits); CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                 "%3d. %s   with eps=%.2e ... %3d KSP (last) iters and %3d Newton iters\n",
                 m+1,SNESConvergedReasons[reason],kspits,its,user.eps);CHKERRQ(ierr);
      if (reason < 0) {
          if ((user.divergetryagain == PETSC_TRUE) && (user.eps > 0.0)) {
              ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "         ... try again w eps *= 1.5\n");CHKERRQ(ierr);
              //ierr = SNESDestroy(&snes); CHKERRQ(ierr);
              //ierr = SNESboot(&snes,&user); CHKERRQ(ierr);
              user.eps = 1.5 * eps_sched[m];
              ierr = VecCopy(H,Htry); CHKERRQ(ierr);
              //ierr = ExplicitStepSmoother(Htry,&user); CHKERRQ(ierr);
              ierr = SNESAttempt(snes,Htry,&its,&reason);CHKERRQ(ierr);
              ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
              ierr = KSPGetIterationNumber(ksp,&kspits); CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "     %s AGAIN  eps=%.2e ... %3d KSP (last) iters and %3d Newton iters\n",
                     SNESConvergedReasons[reason],user.eps,kspits,its);CHKERRQ(ierr);
          }
          if (reason < 0) {
              user.eps = (m>0 ? eps_sched[m-1] : eps_sched[0]); // record last successful eps
              break;
          }
      }
      ierr = VecCopy(Htry,H); CHKERRQ(ierr);
      ierr = StdoutReport(H,&user); CHKERRQ(ierr);
  }
  gettimeofday(&user.endtime, NULL);

  if (user.history == PETSC_TRUE) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
          "writing history.txt to %s ...\n", user.figsprefix); CHKERRQ(ierr);
      ierr = WriteHistoryFile(H,"history.txt",argc,argv,&user); CHKERRQ(ierr);
  }

  if (user.dump == PETSC_TRUE) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
          "writing x.dat,y.dat,b.dat,m.dat,Hexact.dat,H.dat to %s ...\n",
          user.figsprefix); CHKERRQ(ierr);
      ierr = DumpToFiles(H,&user); CHKERRQ(ierr);
  }

  ierr = VecDestroy(&user.m);CHKERRQ(ierr);
  ierr = VecDestroy(&user.b);CHKERRQ(ierr);
  ierr = VecDestroy(&user.Hexact);CHKERRQ(ierr);
  ierr = VecDestroy(&Htry);CHKERRQ(ierr);
  ierr = VecDestroy(&H);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = DMDestroy(&user.quadda);CHKERRQ(ierr);
  ierr = DMDestroy(&user.da);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
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


// averages two gradients; only used in true Mahaffy
Grad gradav(Grad g1, Grad g2) {
  Grad gav;
  gav.x = (g1.x + g2.x) / 2.0;
  gav.y = (g1.y + g2.y) / 2.0;
  return gav;
}


/* For call-back.
Evaluate residual FF on local process patch:
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
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info,PetscScalar **H,PetscScalar **FF,
                                 AppCtx *user) {
  PetscErrorCode  ierr;
  const PetscReal dx = user->dx,
                  dy = dx,
                  upmin = (user->upwindfull == PETSC_TRUE) ? 0.0 : 1.0/4.0,
                  upmax = (user->upwindfull == PETSC_TRUE) ? 1.0 : 3.0/4.0,
                  // lengths & dirs of segments on bdry of control vol:
                  coeff[8] = {dy/2, dx/2, dx/2, -dy/2, -dy/2, -dx/2, -dx/2, dy/2},
                  // local (element-wise) coords of quadrature points for M*
                  locx[4]     = { dx/2.0, 3.0*dx/4.0,     dx/2.0, dx/4.0},
                  locy[4]     = { dy/4.0,     dy/2.0, 3.0*dy/4.0, dy/2.0},
                  locxtrue[4] = { dx/2.0,         dx,     dx/2.0,    0.0},
                  locytrue[4] = {    0.0,     dy/2.0,         dy, dy/2.0},
                  locxnbr[4]  = { dx/2.0,        0.0,     dx/2.0,     dx},
                  locynbr[4]  = {     dy,     dy/2.0,        0.0, dy/2.0};
  const PetscInt  jnbr[4] = { 0,  1,  0, -1},
                  knbr[4] = {-1,  0,  1,  0};
  PetscInt        j, k;
  PetscReal       **am, **ab, ***aq, He[4];
  Grad            gH[4], gb[4];
  Vec             bloc, qloc;

  PetscFunctionBeginUser;
  user->maxD = 0.0;

  // loop over locally-owned elements, including ghosts, forcing or checking
  //     nonnegativity of thickness
  if (user->forceadmissible == PETSC_TRUE) {
      for (k = info->ys-1; k < info->ys + info->ym; k++) {
          for (j = info->xs-1; j < info->xs + info->xm; j++) {
              H[k][j] = PetscAbsReal(H[k][j]);
          }
      }
  } else if (user->checkadmissible == PETSC_TRUE) {
      for (k = info->ys-1; k < info->ys + info->ym; k++) {
          for (j = info->xs-1; j < info->xs + info->xm; j++) {
              if (H[k][j] < 0.0)
                  SETERRQ3(PETSC_COMM_WORLD,1,"ERROR: inadmissible H[%d][%d] = %.3e < 0\n",k,j,H[k][j]);
          }
      }
  }

  // need stencil width on b and locally-computed q
  ierr = DMCreateLocalVector(user->da,&bloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(user->da,user->b,INSERT_VALUES,bloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(user->da,user->b,INSERT_VALUES,bloc);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(user->quadda,&qloc);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(user->da, bloc, &ab);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->m, &am);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(user->quadda, qloc, &aq);CHKERRQ(ierr);
  // loop over locally-owned elements, including ghosts, to get fluxes at
  // c = 0,1,2,3 points in element;  note start at (xs-1,ys-1)
  for (k = info->ys-1; k < info->ys + info->ym; k++) {
      for (j = info->xs-1; j < info->xs + info->xm; j++) {
          PetscInt c;
          // evaluate gradients and (non-upwinded) thicknesses
          if (user->mtrue == PETSC_TRUE) {  // true Mahaffy method
              // this implementation is at least a factor of two inefficient
              Grad gH_nbr[4], gb_nbr[4];
              for (c=0; c<4; c++) {
                  He[c] = fieldatpt(j,k,locxtrue[c],locytrue[c],H, user);
                  gH[c] = gradfatpt(j,k,locxtrue[c],locytrue[c],H, user);
                  gb[c] = gradfatpt(j,k,locxtrue[c],locytrue[c],ab,user);
                  gH_nbr[c] = gradfatpt(j+jnbr[c],k+knbr[c],locxnbr[c],locynbr[c],H, user);
                  gb_nbr[c] = gradfatpt(j+jnbr[c],k+knbr[c],locxnbr[c],locynbr[c],ab,user);
                  gH[c] = gradav(gH[c],gH_nbr[c]);
                  gb[c] = gradav(gb[c],gb_nbr[c]);
              }
          } else {  // M* method, with or without upwind
              for (c=0; c<4; c++) {
                  He[c] = fieldatpt(j,k,locx[c],locy[c],H, user);
                  gH[c] = gradfatpt(j,k,locx[c],locy[c],H, user);
                  gb[c] = gradfatpt(j,k,locx[c],locy[c],ab,user);
              }
          }
          // evaluate fluxes
          if (user->noupwind == PETSC_TRUE) {  // non-upwinding methods
              for (c=0; c<4; c++)
                  aq[k][j][c] = getflux(gH[c],gb[c],He[c],He[c],xdire[c],user);
          } else {  // M* method including first-order upwinding on grad b part
              PetscReal Hup, locxup[4], locyup[4];
              for (c=0; c<4; c++) {  locxup[c] = locx[c];  locyup[c] = locy[c];  }  // copy
              locxup[0] = (gb[0].x <= 0.0) ? upmin*dx : upmax*dx;
              locyup[1] = (gb[1].y <= 0.0) ? upmin*dy : upmax*dy;
              locxup[2] = (gb[2].x <= 0.0) ? upmin*dx : upmax*dx;
              locyup[3] = (gb[3].y <= 0.0) ? upmin*dy : upmax*dy;
              for (c=0; c<4; c++) {
                  Hup = fieldatpt(j,k,locxup[c],locyup[c],H,user);
                  aq[k][j][c] = getflux(gH[c],gb[c],He[c],Hup,xdire[c],user);
              }
          }
      }
  }
  // loop over nodes, not including ghosts, to get residual from quadature over
  // s = 0,1,...,7 points on boundary of control volume (rectangle) around node
  for (k=info->ys; k<info->ys+info->ym; k++) {
      for (j=info->xs; j<info->xs+info->xm; j++) {
          PetscInt s;
          // This is the integral over the control volume boundary using two quadrature
          // points on each side of of the four sides of the rectangular control volume.
          // For M*: two instances of midpoint rule on each side, with two different values of aq[][] per side
          // For true Mahaffy: ditto, but the two values of aq[][] per side are actually the same.
          FF[k][j] = - am[k][j] * dx * dy;
          for (s=0; s<8; s++)
            FF[k][j] += coeff[s] * aq[k+ke[s]][j+je[s]][ce[s]];
      }
  }
  ierr = DMDAVecRestoreArray(user->da, user->m, &am);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, bloc, &ab);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(user->quadda, qloc, &aq);CHKERRQ(ierr);

  ierr = VecDestroy(&bloc); CHKERRQ(ierr);
  ierr = VecDestroy(&qloc); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode ExplicitStepSmoother(Vec H, AppCtx *user) {
    PetscErrorCode  ierr;
    DMDALocalInfo   info;
    PetscReal       **aHloc, **aR, **aH, maxD, tmpeps, deltat, mu;
    PetscInt        j, k;
    Vec             Hloc, R;

    PetscFunctionBeginUser;
    // generate ghosted version of thickness H
    ierr = DMCreateLocalVector(user->da, &Hloc); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->da,H,INSERT_VALUES,Hloc); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da,H,INSERT_VALUES,Hloc); CHKERRQ(ierr);

    // prepare to get residual from FormFunctionLocal()
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(user->da, &R); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, Hloc, &aHloc);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, H, &aH);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, R, &aR);CHKERRQ(ierr);

    // get eps=0 residual corresponding to H
    //user->eps = 0.0;
    tmpeps = user->eps; // set aside current value of eps
    ierr = FormFunctionLocal(&info,aHloc,aR,user);CHKERRQ(ierr); // computes aR and user->maxD
    user->eps = tmpeps; // restore it

    // based on maxD we can take a max-principle stable explicit step
    //dxinvsum = 1.0 / ((1.0 / user->dx) + (1.0 / user->dy));  <- replaces dx/2
    ierr = MPI_Allreduce(&user->maxD,&maxD,1,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD); CHKERRQ(ierr);
    deltat = (user->dx / 2.0) / (2.0 * maxD);

    // take step aH, which is in-place since we have residual
    // mu = deltat / (user->dx * user->dy);
    mu = deltat / (user->dx * user->dx);
    for (k=info.ys; k<info.ys+info.ym; k++) {
        for (j=info.xs; j<info.xs+info.xm; j++) {
            aH[k][j] = PetscMax(0.0, aH[k][j] - mu * aR[k][j]);
        }
    }

    // clean up
    ierr = DMDAVecRestoreArray(user->da, H, &aH);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, R, &aR);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, Hloc, &aHloc);CHKERRQ(ierr);
    ierr = VecDestroy(&Hloc); CHKERRQ(ierr);
    ierr = VecDestroy(&R); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


PetscErrorCode ProcessOptions(AppCtx *user) {
  PetscErrorCode ierr;
  PetscBool      domechosen;
  char           histprefix[512];

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"mah_","options to mahaffy","");CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-bedstep", "use bedrock step exact solution by Jarosh, Schoof, Anslow (2013)",
      NULL,user->bedstep,&user->bedstep,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-checkadmissible", "in FormFunctionLocal(), stop if H < 0.0",
      NULL,user->checkadmissible,&user->checkadmissible,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-dome", "use dome exact solution by Bueler (2003) [default]",
      NULL,user->dome,&user->dome,&domechosen);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-divergetryagain", "on SNES diverge, try again with eps *= 1.5",
      NULL,user->divergetryagain,&user->divergetryagain,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString(
      "-dump", "dump fields into PETSc binary files [x,y,b,m,H,Herror].dat with this prefix",
      NULL,user->figsprefix,user->figsprefix,512,&user->dump); CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-D0", "representative value in m^2/s of diffusivity: D0 ~= D(H,|grad s|)",
      NULL,user->D0,&user->D0,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-forceadmissible", "in FormFunctionLocal(), set H = abs(H)",
      NULL,user->forceadmissible,&user->forceadmissible,NULL);CHKERRQ(ierr);
  strcpy(histprefix,"PREFIX/");
  ierr = PetscOptionsString(
      "-history", "write history.txt file with this prefix (also written under -mah_dump)",
      NULL,histprefix,histprefix,512,&user->history); CHKERRQ(ierr);
  ierr = PetscOptionsInt(
      "-Neps", "levels in schedule of eps regularization/continuation",
      NULL,user->Neps,&user->Neps,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-noupwind", "do not upwind ''W H^{n+2}'' term in flux formula  q = - D grad H + W H^{n+2}",
      NULL,user->noupwind,&user->noupwind,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString(
      "-read", "read grid and data from special-format PETSc binary file; see README.md",
      NULL,user->readname,user->readname,512,&user->read); CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-showdata", "use PETSc X viewers to show b, m, and exact/observed H",
      NULL,user->showdata,&user->showdata,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-swapxy", "swap coordinates x and y when building bedrock step exact solution",
      NULL,user->swapxy,&user->swapxy,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-true", "use true Mahaffy method, not default M*",
      NULL,user->mtrue,&user->mtrue,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-upwindfull", "when upwinding ''W H^{n+2}'' term, go all the way to element edge",
      NULL,user->upwindfull,&user->upwindfull,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  // enforce consistency of cases
  if (user->read == PETSC_TRUE) {
      if (domechosen == PETSC_TRUE) {
        SETERRQ(PETSC_COMM_WORLD,1,"ERROR option conflict: both -mah_dome and -mah_read not allowed\n");
      } else
        user->dome = PETSC_FALSE;    // here user has chosen -mah_read and not -mah_dome
  }
  if (user->bedstep == PETSC_TRUE) {
      if (domechosen == PETSC_TRUE) {
        SETERRQ(PETSC_COMM_WORLD,2,"ERROR option conflict: both -mah_dome and -mah_bedstep not allowed\n");
      } else
        user->dome = PETSC_FALSE;    // here user has chosen -mah_bedstep and not -mah_dome
  }
  if ((user->read == PETSC_TRUE) && (user->bedstep == PETSC_TRUE)) {
      SETERRQ(PETSC_COMM_WORLD,3,"ERROR option conflict: both -mah_bedstep and -mah_read not allowed\n");
  }
  if (user->mtrue == PETSC_TRUE) {
      user->noupwind = PETSC_TRUE;
  }
  if (user->dump == PETSC_TRUE) {
      user->history = PETSC_TRUE;
  } else if (user->history == PETSC_TRUE) {
      strcpy(user->figsprefix,histprefix);
  }
  PetscFunctionReturn(0);
}


// start a SNESVI
PetscErrorCode SNESboot(SNES *s, AppCtx* user) {
  PetscErrorCode ierr;
  ierr = SNESCreate(PETSC_COMM_WORLD,s);CHKERRQ(ierr);
  ierr = SNESSetDM(*s,user->da);CHKERRQ(ierr);
  ierr = DMDASNESSetFunctionLocal(user->da,INSERT_VALUES,
              (DMDASNESFunction)FormFunctionLocal,user); CHKERRQ(ierr);
  ierr = SNESVISetComputeVariableBounds(*s,&FormBounds);CHKERRQ(ierr);
  ierr = SNESSetType(*s,SNESVINEWTONRSLS);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(*s);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


// try a SNES solve; H is modified
PetscErrorCode SNESAttempt(SNES snes, Vec H, PetscInt *its, SNESConvergedReason *reason) {
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = SNESSolve(snes, NULL, H); CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(snes,its);CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(snes,reason);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

