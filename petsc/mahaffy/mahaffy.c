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
"Compares Mahaffy* (default) and true Mahaffy schemes.\n"
"Uses SNESVI (-snes_type vinewtonrsls) with constraint  H(x,y) >= 0.\n\n"
"Two cases: (1) flat bed case where analytical solution is known,\n"
"           (2) -mah_read foo.dat reads b and m from PETSc binary file.\n\n"
"For case (1) usage read usage comment at start of mahaffy.c.\n"
"For case (2) usage see README.md.\n\n";

/* Usage help:
   ./mahaffy -help |grep mah_

Use one of these three:
   ./mahaffy -mah_dome         # default
   ./mahaffy -mah_bedstep
   ./mahaffy -mah_read foo.dat # see README.md for Greenland example

Solution options:
   ./mahaffy -mah_true         # use true Mahaffy
   ./mahaffy -mah_upwind       # a kind of first-order upwinding on  grad b  part of flux
   ./mahaffy -mah_Neps 4       # don't go all the way on continuation
   ./mahaffy -mah_D0 10.0      # change constant diffusivity, used in continuation, to other value than
                               # default = 1.0;  big may be good for convergence, esp. w upwinding
   ./mahaffy -mah_divergetryagain # on SNES diverge, try again with eps *= 1.5

   ./mahaffy -da_refine 1 -snes_vi_monitor  # widen screen to see SNESVI monitor output

Successes:
   ./mahaffy -mah_bedstep -mah_D0 0.01 -da_refine 4 -mah_upwind -snes_max_it 1000  # needs 604 SNES iterations!
   mpiexec -n 6 ./mahaffy -pc_type mg -da_refine 5 -snes_max_it 200 -snes_monitor


Divergence and error cases:
   ./mahaffy -da_refine 2      # divergence at 9
   ./mahaffy -da_refine 3      # divergence at 9

Compare methods:
   mpiexec -n 6 ./mahaffy -da_refine 4 -mah_Neps 10
   mpiexec -n 6 ./mahaffy -da_refine 4 -mah_Neps 10 -mah_true   # line 505 of virs.c

High-res success (also -da_refine 7 works):
   mpiexec -n 6 ./mahaffy -da_refine 6 -pc_type asm -sub_pc_type lu -snes_max_it 200

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


extern PetscErrorCode FormBounds(SNES,Vec,Vec);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*,PetscScalar**,PetscScalar**,AppCtx*);
// extern PetscErrorCode FormJacobianLocal(DMDALocalInfo*,PetscScalar*,Mat,Mat,AppCtx*); //see layer.c

extern PetscErrorCode SNESboot(SNES*,AppCtx*);
extern PetscErrorCode SNESAttempt(SNES,Vec,PetscInt*,SNESConvergedReason*);

extern PetscErrorCode ProcessOptions(AppCtx*);
extern PetscErrorCode StateReport(Vec,DMDALocalInfo*,AppCtx*);


int main(int argc,char **argv) {
  PetscErrorCode      ierr;
  SNES                snes;
  Vec                 H, Htry;
  AppCtx              user;
  DMDALocalInfo       info;

  PetscInitialize(&argc,&argv,(char*)0,help);

  // default settings of parameters
  user.Nx     = -12;        // so DMDACreate2d() defaults to 12x12
  user.Ny     = -12;

  user.n      = 3.0;
  user.g      = 9.81;       // m/s^2
  user.rho    = 910.0;      // kg/m^3
  user.secpera= 31556926.0;
  user.A      = 1.0e-16/user.secpera; // = 3.17e-24  1/(Pa^3 s); EISMINT I value
  user.Gamma  = 2.0 * PetscPowReal(user.rho*user.g,user.n) * user.A / (user.n+2.0);

  user.eps    = 0.0;
  user.D0     = 1.0;        // m^2 / s
  user.Neps   = 13;
  user.mtrue  = PETSC_FALSE;
  user.upwind = PETSC_FALSE;

  user.read   = PETSC_FALSE;
  user.dome   = PETSC_TRUE;  // defaults to this case
  user.bedstep= PETSC_FALSE;
  user.swapxy = PETSC_FALSE;
  user.divergetryagain = PETSC_FALSE;

  user.showdata= PETSC_FALSE;
  user.dump   = PETSC_FALSE;
  strcpy(user.figsprefix,"PREFIX/");  // dummies improve "mahaffy -help" appearance
  strcpy(user.readname,"FILENAME");

  ierr = ProcessOptions(&user); CHKERRQ(ierr);

  if (user.read == PETSC_TRUE) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
          "reading dimensions and grid spacing from %s ...\n",user.readname); CHKERRQ(ierr);
      ierr = ReadDimensions(&user); CHKERRQ(ierr);  // sets user.[Nx,Ny,Lx,Ly,dx]
  } else {
      if (user.dome == PETSC_TRUE) {
          user.Lx     = 900.0e3;    // m
          user.Ly     = 900.0e3;    // m
      } else {
          user.Lx     = 30.0e3;    // m
          user.Ly     = 30.0e3;    // m
      }
  }

  // this DMDA is used for scalar fields on nodes; cell-centered grid
  ierr = DMDACreate2d(PETSC_COMM_WORLD,
                      DM_BOUNDARY_PERIODIC,DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX,
                      user.Nx,user.Ny,PETSC_DECIDE,PETSC_DECIDE,  // default grid if Nx<0, Ny<0
                      1, 1,                                       // dof=1,  stencilwidth=1
                      NULL,NULL,&user.da);
  ierr = DMSetApplicationContext(user.da, &user);CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(user.da,&info); CHKERRQ(ierr);
  if (user.read == PETSC_FALSE) {
      user.dx = 2.0 * user.Lx / (PetscReal)(info.mx);
  }
  ierr = DMDASetUniformCoordinates(user.da, -user.Lx+user.dx/2, user.Lx+user.dx/2,
                                            -user.Ly+user.dx/2, user.Ly+user.dx/2,
                                   0.0,1.0); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      "solving on [-Lx,Lx]x[-Ly,Ly] with Lx=%.3f km and Ly=%.3f km,  %d x %d grid,  spacing dx = %.6f km ...\n",
      user.Lx/1000.0,user.Ly/1000.0,info.mx,info.my,user.dx/1000.0); CHKERRQ(ierr);

  // this DMDA is used for evaluating flux components at quad points on elements
  ierr = DMDACreate2d(PETSC_COMM_WORLD,
                      DM_BOUNDARY_PERIODIC,DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX,
                      info.mx,info.my,PETSC_DECIDE,PETSC_DECIDE,
                      4, 1,                               // dof=4,  stencilwidth=1
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
  ierr = VecScale(H,1500.0*user.secpera); CHKERRQ(ierr);  // FIXME make user.initializemagic
  // alternatives:
  //ierr = [VERIF]ExactThickness(H,&user);CHKERRQ(ierr);
  //ierr = VecSet(H,0.0); CHKERRQ(ierr);

  if (user.mtrue == PETSC_TRUE) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,    "solving by true Mahaffy method ...\n"); CHKERRQ(ierr);
  } else {
      if (user.upwind == PETSC_TRUE) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"solving by Mahaffy*-with-upwinding (Mahaffy++) method ...\n"); CHKERRQ(ierr);
      } else {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"solving by Mahaffy* method (without upwinding) ...\n"); CHKERRQ(ierr);
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
      ierr = StdoutReport(H,&info,&user); CHKERRQ(ierr);
  }
  gettimeofday(&user.endtime, NULL);

  if (user.dump == PETSC_TRUE) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
          "writing history.txt,x.dat,y.dat,b.dat,m.dat,Hexact.dat,H.dat to %s ...\n",
          user.figsprefix); CHKERRQ(ierr);
      ierr = WriteHistoryFile(H,"history.txt",argc,argv,&info,&user); CHKERRQ(ierr);
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


typedef struct {
    PetscReal x,y;
} Grad;

/* the first formula for SIA: compute delta, a factor of diffusivity D and
pseudo-velocity W, from values of thickness and surface gradient
note:
   delta = Gamma |grad s|^{n-1}
where s = H+b; also applies power-regularization part of continuation scheme */
PetscReal getdelta(const Grad gH, const Grad gb, const AppCtx *user) {
    const PetscReal eps = user->eps,
                    n   = (1.0 - eps) * user->n + eps * 1.0,
                    sx  = gH.x + gb.x,
                    sy  = gH.y + gb.y;
    return user->Gamma * PetscPowReal(sx*sx + sy*sy,(n-1.0)/2);
}


enum Direction  {X, Y};
typedef enum Direction Direction;

/* the second formula for SIA: evaluate the flux from gradients and thickness;
note:
   D = delta * H^{n+2}      (positive scalar)
   W = - delta * grad b     (vector)
so
   q = - D grad H + W H^{n+2}
(and also q = - D grad s) where  H = H_{j*,k*}  at point (x_{j*}, y_{k*});
note  W = (W.x,W.y)  and  q = (q.x,q.y)  in code;
in this upwinding case, for  j <= j* <= j+1,  k <= k* <= k+1,
   q.x = - D gH.x + W.x^+ H_{j,k*}^{n+2} + W.x^- H_{j+1,k*}^{n+2}
   q.y = - D gH.y + W.y^+ H_{j*,k}^{n+2} + W.y^- H_{j*,k+1}^{n+2}
and so if  W.x >= 0  and  W.y >= 0  the input should have inputs
   H = H_{j*,k*},  Hxup = H_{j,k*},  Hyup = H_{j*,k},
this method also applies diffusivity- and power-regularization parts of the
continuation scheme, and it updates maximum of diffusivity */
PetscReal getfluxUP(const Grad gH, const Grad gb,
                    const PetscReal H, const PetscReal Hup,
                    const Direction dir,
                    AppCtx *user) {
  const PetscReal eps   = user->eps,
                  n     = (1.0 - eps) * user->n + eps * 1.0,
                  delta = getdelta(gH,gb,user),
                  D     = (1.0-eps) * delta * PetscPowReal(H,n+2.0) + eps * user->D0;
  user->maxD = PetscMax(user->maxD, D);
  if (dir == X) {
      const PetscReal  Wx = - delta * gb.x;
      return - D * gH.x + Wx * PetscPowReal(Hup,n+2.0);
  } else {
      const PetscReal  Wy = - delta * gb.y;
      return - D * gH.y + Wy * PetscPowReal(Hup,n+2.0);
  }
}


// the non-upwinding form of the method, i.e. pure Mahaffy*
PetscReal getflux(const Grad gH, const Grad gb, // compute with gradsatpt() first
                  const PetscReal H,            // compute with thicknessatpt() first
                  const Direction dir,
                  AppCtx *user) {
  return getfluxUP(gH,gb,H,H,dir,user);
}


/* first of two operations with Q1 interpolants on an element:
   evaluate gradients of the thickness and bed elevation at any point (x,y) on element, using corner values of H and b */
PetscErrorCode gradsatpt(PetscInt j, PetscInt k,         // (j,k) is the element (by lower-left corner)
                         PetscReal locx, PetscReal locy, // = (x,y) coords in element
                         PetscReal **H, PetscReal **b,   // H[k][j] and b[k][j] are node values
                         const AppCtx *user,
                         Grad *gradH, Grad *gradb) {

  const PetscReal dx = user->dx, dy = dx,
                  x[4]  = {1.0 - locx / dx, locx / dx,       1.0 - locx / dx, locx / dx},
                  gx[4] = {- 1.0 / dx,      1.0 / dx,        - 1.0 / dx,      1.0 / dx},
                  y[4]  = {1.0 - locy / dy, 1.0 - locy / dy, locy / dy,       locy / dy},
                  gy[4] = {- 1.0 / dy,      - 1.0 / dy,      1.0 / dy,        1.0 / dy};
  if ((H == NULL) || (b == NULL)) {
      SETERRQ(PETSC_COMM_WORLD,1,"ERROR: illegal NULL ptr in gradsatpt() ...\n");
  }
  gradH->x = gx[0] * y[0] * H[k][j]   + gx[1] * y[1] * H[k][j+1]
           + gx[2] * y[2] * H[k+1][j] + gx[3] * y[3] * H[k+1][j+1];
  gradH->y =  x[0] *gy[0] * H[k][j]   +  x[1] *gy[1] * H[k][j+1]
           +  x[2] *gy[2] * H[k+1][j] +  x[3] *gy[3] * H[k+1][j+1];
  gradb->x = gx[0] * y[0] * b[k][j]   + gx[1] * y[1] * b[k][j+1]
           + gx[2] * y[2] * b[k+1][j] + gx[3] * y[3] * b[k+1][j+1];
  gradb->y =  x[0] *gy[0] * b[k][j]   +  x[1] *gy[1] * b[k][j+1]
           +  x[2] *gy[2] * b[k+1][j] +  x[3] *gy[3] * b[k+1][j+1];
  PetscFunctionReturn(0);
}


/* second of two operations with Q1 interpolants on an element:
   evaluate the thickness H at any point (x,y) on element, using corner values of H */
PetscErrorCode thicknessatpt(PetscInt j, PetscInt k,         // (j,k) is the element (by lower-left corner)
                             PetscReal locx, PetscReal locy, // = (x,y) coords in element
                             PetscReal **H,                  // H[k][j] are node values
                             AppCtx *user, PetscReal *Hatpt) {
  const PetscReal dx = user->dx,  dy = dx,
                  x[4]  = {1.0 - locx / dx, locx / dx,       1.0 - locx / dx, locx / dx},
                  y[4]  = {1.0 - locy / dy, 1.0 - locy / dy, locy / dy,       locy / dy};
  if (H == NULL) {
      SETERRQ(PETSC_COMM_WORLD,1,"ERROR: illegal NULL ptr in thicknessatpt() ...\n");
  }
  *Hatpt = x[0] * y[0] * H[k][j] + x[1] * y[1] * H[k][j+1] + x[2] * y[2] * H[k+1][j] + x[3] * y[3] * H[k+1][j+1];
  PetscFunctionReturn(0);
}


// averages two gradients; only used in true Mahaffy
Grad gradav(Grad g1, Grad g2) {
  Grad gav;
  gav.x = (g1.x + g2.x) / 2.0;
  gav.y = (g1.y + g2.y) / 2.0;
  return gav;
}


/* on element shown, indexed by (j,k) node at lower-left corner @,
 FluxQuad holds x-components at * and y-components at %
   ---------------
  |       |       |
  |       *xn     |
  |  yw   |  ye   |
  |---%--- ---%---|
  |       |       |
  |       *xs     |
  |       |       |
  @---------------
(j,k)
*/
typedef struct {
    PetscReal xn, xs, ye, yw;
} FluxQuad;


/* for call-back: evaluate residual FF(x) on local process patch */
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info,PetscScalar **H,PetscScalar **FF,
                                 AppCtx *user) {
  PetscErrorCode  ierr;
  const PetscReal dx = user->dx, dy = dx;
  PetscInt        j, k;
  PetscReal       **am, **ab;
  FluxQuad        **aq;
  Grad            gHn, gHs, gHe, gHw,
                  gbn, gbs, gbe, gbw;
  PetscReal       Hn, Hs, He, Hw;
  PetscReal       qnx, qsx, qey, qwy;
  Vec             bloc, qloc;

  PetscFunctionBeginUser;
  user->maxD = 0.0;

  ierr = DMCreateLocalVector(user->da,&bloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(user->da,user->b,INSERT_VALUES,bloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(user->da,user->b,INSERT_VALUES,bloc);CHKERRQ(ierr);

  ierr = DMCreateLocalVector(user->quadda,&qloc);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(user->da, bloc, &ab);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->m, &am);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->quadda, qloc, &aq);CHKERRQ(ierr);
  // loop over locally-owned elements, including ghosts, to get fluxes
  // note start at (xs-1,ys-1)
  for (k = info->ys-1; k < info->ys + info->ym; k++) {
      for (j = info->xs-1; j < info->xs + info->xm; j++) {
          // get gradients and (non-upwinded) thicknesses
          if (user->mtrue == PETSC_FALSE) {  // Mahaffy* or Mahaffy++ method
              // above-center point in element
              ierr =     gradsatpt(j,k,     dx/2.0, 3.0*dy/4.0,H,ab,user,&gHn,&gbn); CHKERRQ(ierr);
              ierr = thicknessatpt(j,k,     dx/2.0, 3.0*dy/4.0,H,   user,&Hn); CHKERRQ(ierr);
              // below-center point in element
              ierr =     gradsatpt(j,k,     dx/2.0,     dy/4.0,H,ab,user,&gHs,&gbs); CHKERRQ(ierr);
              ierr = thicknessatpt(j,k,     dx/2.0,     dy/4.0,H,   user,&Hs); CHKERRQ(ierr);
              // right-of-center point in element
              ierr =     gradsatpt(j,k, 3.0*dx/4.0,     dy/2.0,H,ab,user,&gHe,&gbe); CHKERRQ(ierr);
              ierr = thicknessatpt(j,k, 3.0*dx/4.0,     dy/2.0,H,   user,&He); CHKERRQ(ierr);
              // left-of-center point in element
              ierr =     gradsatpt(j,k,     dx/4.0,     dy/2.0,H,ab,user,&gHw,&gbw); CHKERRQ(ierr);
              ierr = thicknessatpt(j,k,     dx/4.0,     dy/2.0,H,   user,&Hw ); CHKERRQ(ierr);
          } else {  // true Mahaffy method
                    // this implementation is at least a factor of two inefficient
              // center-top point in element
              ierr =     gradsatpt(j,k,   dx/2.0, dy,H,ab,user,&gHn,&gbn); CHKERRQ(ierr);
              if (k < info->ys + info->ym - 1) {
                Grad gHn_nbr, gbn_nbr;
                ierr = gradsatpt(j,k+1, dx/2.0, 0.0,H,ab,user,&gHn_nbr,&gbn_nbr); CHKERRQ(ierr);
                gHn  = gradav(gHn,gHn_nbr);
                gbn  = gradav(gbn,gbn_nbr);
              }
              ierr = thicknessatpt(j,k,   dx/2.0, dy,H,   user,&Hn); CHKERRQ(ierr);
              // center-bottom point in element
              ierr =     gradsatpt(j,k,   dx/2.0, 0.0,H,ab,user,&gHs,&gbs); CHKERRQ(ierr);
              if (k > info->ys - 1) {
                Grad gHs_nbr, gbs_nbr;
                ierr = gradsatpt(j,k-1, dx/2.0, dy,H,ab,user,&gHs_nbr,&gbs_nbr); CHKERRQ(ierr);
                gHs  = gradav(gHs,gHs_nbr);
                gbs  = gradav(gbs,gbs_nbr);
              }
              ierr = thicknessatpt(j,k,   dx/2.0, 0.0,H,   user,&Hs); CHKERRQ(ierr);
              // center-right point in element
              ierr =     gradsatpt(j,k,   dx,  dy/2.0,H,ab,user,&gHe,&gbe); CHKERRQ(ierr);
              if (j < info->xs + info->xm - 1) {
                Grad gHe_nbr, gbe_nbr;
                ierr = gradsatpt(j+1,k, 0.0, dy/2.0,H,ab,user,&gHe_nbr,&gbe_nbr); CHKERRQ(ierr);
                gHe  = gradav(gHe,gHe_nbr);
                gbe  = gradav(gbe,gbe_nbr);
              }
              ierr = thicknessatpt(j,k,   dx,  dy/2.0,H,   user,&He); CHKERRQ(ierr);
              // center-left point in element
              ierr =     gradsatpt(j,k,   0.0, dy/2.0,H,ab,user,&gHw,&gbw); CHKERRQ(ierr);
              if (j > info->xs - 1) {
                Grad gHw_nbr, gbw_nbr;
                ierr = gradsatpt(j-1,k, dx,  dy/2.0,H,ab,user,&gHw_nbr,&gbw_nbr); CHKERRQ(ierr);
                gHw  = gradav(gHw,gHw_nbr);
                gbw  = gradav(gbw,gbw_nbr);
              }
              ierr = thicknessatpt(j,k,   0.0, dy/2.0,H,   user,&Hw); CHKERRQ(ierr);
          }
          // evaluate fluxes
          if (user->upwind == PETSC_TRUE) {
              // Mahaffy++ method = Mahaffy* + first-order upwinding on grad b part
              PetscReal Hnup, Hsup, Heup, Hwup;
              if (gbn.x <= 0.0) {  // W.x >= 0 case
                  ierr = thicknessatpt(j,k,    dx/4.0,3.0*dy/4.0,H,user,&Hnup); CHKERRQ(ierr);
              } else {
                  ierr = thicknessatpt(j,k,3.0*dx/4.0,3.0*dy/4.0,H,user,&Hnup); CHKERRQ(ierr);
              }
              if (gbs.x <= 0.0) {  // W.x >= 0 case
                  ierr = thicknessatpt(j,k,    dx/4.0,    dy/4.0,H,user,&Hsup); CHKERRQ(ierr);
              } else {
                  ierr = thicknessatpt(j,k,3.0*dx/4.0,    dy/4.0,H,user,&Hsup); CHKERRQ(ierr);
              }
              if (gbe.y <= 0.0) {  // W.y >= 0 case
                  ierr = thicknessatpt(j,k,3.0*dx/4.0,    dy/4.0,H,user,&Heup); CHKERRQ(ierr);
              } else {
                  ierr = thicknessatpt(j,k,3.0*dx/4.0,3.0*dy/4.0,H,user,&Heup); CHKERRQ(ierr);
              }
              if (gbw.y <= 0.0) {  // W.y >= 0 case
                  ierr = thicknessatpt(j,k,    dx/4.0,    dy/4.0,H,user,&Hwup); CHKERRQ(ierr);
              } else {
                  ierr = thicknessatpt(j,k,    dx/4.0,3.0*dy/4.0,H,user,&Hwup); CHKERRQ(ierr);
              }
              qnx = getfluxUP(gHn,gbn,Hn,Hnup,X,user);
              qsx = getfluxUP(gHs,gbs,Hs,Hsup,X,user);
              qey = getfluxUP(gHe,gbe,He,Heup,Y,user);
              qwy = getfluxUP(gHw,gbw,Hw,Hwup,Y,user);
          } else { // both true Mahaffy and Mahaffy* methods
              qnx = getflux(gHn,gbn,Hn,X,user);
              qsx = getflux(gHs,gbs,Hs,X,user);
              qey = getflux(gHe,gbe,He,Y,user);
              qwy = getflux(gHw,gbw,Hw,Y,user);
          }
          aq[k][j] = (FluxQuad){qnx,qsx,qey,qwy};
      }
  }
  // loop over nodes to get residual
  // This is the integral over boundary of control volume using two quadrature
  // points on each side of control volume.  For Mahaffy* it is two instances of
  // midpoint rule on each side.  For Mahaffy it is just one, but with averaged
  // gradient calculation at the midpoint of the side.
  for (k=info->ys; k<info->ys+info->ym; k++) {
      for (j=info->xs; j<info->xs+info->xm; j++) {
          FF[k][j] =  (aq[k][j].xs     + aq[k-1][j].xn  ) * dy/2.0
                    + (aq[k][j-1].ye   + aq[k][j].yw    ) * dx/2.0
                    - (aq[k][j-1].xs   + aq[k-1][j-1].xn) * dy/2.0
                    - (aq[k-1][j-1].ye + aq[k-1][j].yw  ) * dx/2.0
                    - am[k][j] * dx * dy;
      }
  }
  ierr = DMDAVecRestoreArray(user->da, user->m, &am);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, bloc, &ab);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->quadda, qloc, &aq);CHKERRQ(ierr);

  ierr = VecDestroy(&bloc); CHKERRQ(ierr);
  ierr = VecDestroy(&qloc); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode ProcessOptions(AppCtx *user) {
  PetscErrorCode ierr;
  PetscBool      domechosen;

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"mah_","options to mahaffy","");CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-bedstep", "use bedrock step exact solution by Jarosh, Schoof, Anslow (2013)",
      NULL,user->bedstep,&user->bedstep,NULL);CHKERRQ(ierr);
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
  ierr = PetscOptionsInt(
      "-Neps", "levels in schedule of eps regularization/continuation",
      NULL,user->Neps,&user->Neps,NULL);CHKERRQ(ierr);
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
      "-true", "use true Mahaffy method, not default Mahaffy*",
      NULL,user->mtrue,&user->mtrue,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-upwind", "upwind ''W H^{n+2}'' term in flux formula  q = - D grad H + W H^{n+2}",
      NULL,user->upwind,&user->upwind,NULL);CHKERRQ(ierr);
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


