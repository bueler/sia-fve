static const char help[] =
"Solves steady ice sheet problem in 2d:\n"
"    div (q^x,q^y) = m,\n"
"    (q^x,q^y) = - Gamma H^{n+2} |grad s|^{n-1} grad s,\n"
"where  H(x,y)  is the ice thickness,  b(x,y)  is the bed elevation,\n"
"and  s(x,y) = H(x,y) + b(x,y)  is surface elevation.\n"
"Note  n > 1  and  Gamma = 2 A (rho g)^n / (n+2).\n"
"Domain is  -L < x < L,  -L < y < L,  with periodic boundary conditions.\n"
"Computed in flat bed case where analytical solution is known.\n"
"Computed by Q1 FVE method with FD evaluation of Jacobian (i.e. no analytical yet).\n"
"Compares Mahaffy* (default) and true Mahaffy schemes.\n"
"Uses SNESVI = -snes_type vinewtonrsls.\n\n";

//   ./mahaffy -help |grep mah_

//   ./mahaffy

//   ./mahaffy -mah_true       // use true Mahaffy
//   ./mahaffy -mah_Neps 4     // don't go all the way on continuation

// widen screen to see SNESVI monitor output:
//   ./mahaffy -da_refine 1 -snes_vi_monitor

// error from:
//   ./mahaffy -da_refine 2  // divergence, but no crash
//   ./mahaffy -da_refine 3  // crash at SNESSolve_VINEWTONRSLS() line 505
                             //     in /home/ed/petsc-maint/src/snes/impls/vi/rs/virs.c
//   ./mahaffy -da_refine 1 -pc_type mg  // error at line 506 of virs.c

// compare methods
//   mpiexec -n 6 ./mahaffy -da_refine 4 -mah_Neps 10
//   mpiexec -n 6 ./mahaffy -da_refine 4 -mah_Neps 10 -mah_true   // again, line 505 of virs.c

// high-res success (also -da_refine 7 works):
//   mpiexec -n 6 ./mahaffy -da_refine 6 -pc_type asm -sub_pc_type lu -snes_max_it 200

// generate .png figs:
//   mkdir foo/
//   ./mahaffy -mah_dump foo/
//   cd foo/
//   python ../figsmahaffy.py

#include <math.h>
#include <petscdmda.h>
#include <petscsnes.h>

typedef struct {
  DM        da, quadda;
  Vec       b,      // the bed elevation
            m,      // the (steady) surface mass balance
            topgread,cmbread;  // only valid if user.read is TRUE and ReadFromBinary() has happened
  PetscReal dx,     // fixed grid spacing; dx = dy
            Lx,     // domain is [-Lx,Lx] x [-Ly,Ly]
            Ly,
            n,      // Glen exponent for SIA flux term
            g,      // acceleration of gravity
            rho,    // ice density
            secpera,// number of seconds in a year
            A,      // ice softness
            Gamma,  // coefficient for SIA flux term
            eps,    // dimensionless regularization
            D0,     // representative value of diffusivity (in regularization)
            exactL, // radius of exact ice sheet
            exactH0;// center thickness of exact ice sheet
  PetscInt  Nx,     // grid has Nx x Ny nodes
            Ny,
            Neps;   // number of levels in regularization/continuation
  PetscBool mtrue,  // use true Mahaffy method instead of Mahaffy* (default)
            read,   // read grid and data from PETSc binary file
            showdata,// show b and m with X viewer
            dump;   // dump fields into ASCII VTK files
  char      figsprefix[512],
            readname[512];
} AppCtx;


extern PetscErrorCode FormBounds(SNES,Vec,Vec);
extern PetscErrorCode SetToExactThickness(Vec,const AppCtx*);
extern PetscErrorCode SetToSMB(Vec,PetscBool,const AppCtx*);
extern PetscErrorCode ProcessOptions(AppCtx*);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*,PetscScalar**,PetscScalar**,AppCtx*);
// extern PetscErrorCode FormJacobianLocal(DMDALocalInfo*,PetscScalar*,Mat,Mat,AppCtx*); //see layer.c
extern PetscErrorCode SNESAttempt(SNES,Vec,PetscInt*,SNESConvergedReason*);
extern PetscErrorCode ErrorReport(Vec,DMDALocalInfo*,AppCtx*);
extern PetscErrorCode ViewToVTKASCII(Vec,const char[],const char[]);
extern PetscErrorCode DumpToFiles(Vec,AppCtx*);
extern PetscErrorCode SNESboot(SNES*,AppCtx*);
extern PetscErrorCode ReadFromBinary(AppCtx*);
extern PetscErrorCode ReshapeReadVecs(AppCtx*);


int main(int argc,char **argv) {
  PetscErrorCode      ierr;
  SNES                snes;
  Vec                 H, Htry;
  AppCtx              user;
  DMDALocalInfo       info;

  PetscInitialize(&argc,&argv,(char*)0,help);

  user.Nx     = -12;        // so DMDACreate2d() defaults to 12x12
  user.Ny     = -12;
  user.Lx     = 900.0e3;    // m
  user.Ly     = 900.0e3;    // m

  user.n      = 3.0;
  user.g      = 9.81;       // m/s^2
  user.rho    = 910.0;      // kg/m^3
  user.secpera= 31556926.0;
  user.A      = 1.0e-16/user.secpera; // = 3.17e-24  1/(Pa^3 s); EISMINT I value
  user.Gamma  = 2.0 * PetscPowReal(user.rho*user.g,user.n) * user.A / (user.n+2.0);

  user.exactL = 750.0e3;    // m
  user.exactH0= 3600.0;     // m

  user.eps    = 0.0;
  user.D0     = 1.0;        // m^2 / s
  user.Neps   = 13;
  user.mtrue  = PETSC_FALSE;

  user.showdata= PETSC_FALSE;
  user.dump   = PETSC_FALSE;
  user.read   = PETSC_FALSE;
  strcpy(user.figsprefix,"PREFIX/");  // dummies improve "mahaffy -help" appearance
  strcpy(user.readname,"FILENAME");

  ierr = ProcessOptions(&user); CHKERRQ(ierr);

  if (user.read == PETSC_TRUE) {
      // sets user.Nx, user.Ny, user.Lx, user.Ly, user.dx and reads serial Vecs topgread,cmbread
      ierr = ReadFromBinary(&user); CHKERRQ(ierr);
  }

  // this DMDA is used for scalar fields on nodes; cell-centered grid
  ierr = DMDACreate2d(PETSC_COMM_WORLD,
                      DM_BOUNDARY_PERIODIC,DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX,
                      user.Nx,user.Ny,PETSC_DECIDE,PETSC_DECIDE,  // default grid if Nx<0, Ny<0
                      1, 1,                               // dof=1,  stencilwidth=1
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
      "solving on [-Lx,Lx]x[-Ly,Ly] with Lx=%.3f km and Ly=%.3f,  %d x %d grid,  spacing dx = %.6f km ...\n",
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

  if (user.read == PETSC_TRUE) {
      ierr = ReshapeReadVecs(&user); CHKERRQ(ierr);
  } else {
      ierr = VecSet(user.b,0.0); CHKERRQ(ierr);
      ierr = SetToSMB(user.m,PETSC_FALSE,&user); CHKERRQ(ierr);
  }

  if (user.showdata == PETSC_TRUE) {
      PetscViewer graphical;
      ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,"bed topography (m)",
                                 PETSC_DECIDE,PETSC_DECIDE,
                                 info.mx,info.my,&graphical); CHKERRQ(ierr);
      ierr = VecView(user.b,graphical); CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&graphical); CHKERRQ(ierr);
      ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,"climatic mass balance (m s-1)",
                                 PETSC_DECIDE,PETSC_DECIDE,
                                 info.mx,info.my,&graphical); CHKERRQ(ierr);
      ierr = VecView(user.m,graphical); CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&graphical); CHKERRQ(ierr);
  }

  // initialize by chop & scale SMB
  ierr = SetToSMB(H,PETSC_TRUE,&user); CHKERRQ(ierr);
  ierr = VecScale(H,1500.0*user.secpera); CHKERRQ(ierr);  // FIXME make user.initializemagic
  // alternatives:
  //ierr = SetToExactThickness(H,&user);CHKERRQ(ierr);
  //ierr = VecSet(H,0.0); CHKERRQ(ierr);

  ierr = SNESboot(&snes,&user); CHKERRQ(ierr);

  KSP        ksp;
  PetscInt   its, kspits, m;
  PetscReal  eps_sched[13] = {1.0,    0.5,    0.2,     0.1,     0.05,
                                       0.02,   0.01,   0.005,   0.002,   0.001,
                                       0.0005, 0.0002, 0.0000};
  SNESConvergedReason reason;
  for (m = 0; m<user.Neps; m++) {
      user.eps = eps_sched[m];
      ierr = VecCopy(H,Htry); CHKERRQ(ierr);
      ierr = SNESAttempt(snes,Htry,&its,&reason);CHKERRQ(ierr);
      ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
      ierr = KSPGetIterationNumber(ksp,&kspits); CHKERRQ(ierr);
      if (reason < 0) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "%3d. DIVERGED   with eps=%.2e ... %3d KSP (last) iters and %3d Newton iters (%s)\n",
                     m+1,kspits,its,user.eps,SNESConvergedReasons[reason]);CHKERRQ(ierr);
          if (user.eps > 0.0) {
              ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "         ... try again w eps *= 2\n");CHKERRQ(ierr);
              //ierr = SNESDestroy(&snes); CHKERRQ(ierr);
              //ierr = SNESboot(&snes,&user); CHKERRQ(ierr);
              user.eps = 2.0 * eps_sched[m];
              ierr = VecCopy(H,Htry); CHKERRQ(ierr);
              ierr = SNESAttempt(snes,Htry,&its,&reason);CHKERRQ(ierr);
              ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
              ierr = KSPGetIterationNumber(ksp,&kspits); CHKERRQ(ierr);
              if (reason < 0) {
                  ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "     DIVERGED AGAIN  eps=%.2e ... %3d KSP (last) iters and %3d Newton iters (%s)\n",
                         user.eps,kspits,its,SNESConvergedReasons[reason]);CHKERRQ(ierr);
                  break;
              }
          }
      }
      if (reason >= 0) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "%3d. converged  with eps=%.2e ... %3d KSP (last) iters and %3d Newton iters (%s)\n",
                     m+1,kspits,its,user.eps,SNESConvergedReasons[reason]);CHKERRQ(ierr);
          ierr = VecCopy(Htry,H); CHKERRQ(ierr);
          if (user.read == PETSC_FALSE) {
              ierr = ErrorReport(H,&info,&user); CHKERRQ(ierr);
          }
      }
  }

  if (user.dump == PETSC_TRUE) {
      ierr = DumpToFiles(H,&user); CHKERRQ(ierr);
  }

  if (user.read == PETSC_TRUE) {
      ierr = VecDestroy(&user.topgread);CHKERRQ(ierr);
      ierr = VecDestroy(&user.cmbread);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&user.m);CHKERRQ(ierr);
  ierr = VecDestroy(&user.b);CHKERRQ(ierr);
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


PetscReal radialcoord(const DMDACoor2d c) {
  PetscReal r;
  r = PetscSqrtReal(c.x * c.x + c.y * c.y);
  if (r < 0.01)
      r = 0.01;  // avoid r=0 singularity
  return r;
}


PetscErrorCode SetToExactThickness(Vec H, const AppCtx *user) {
  PetscErrorCode ierr;

  const PetscReal L  = user->exactL,
                  H0 = user->exactH0,
                  n  = user->n,
                  mm = 1.0 + 1.0 / n,
                  qq = n / (2.0 * n + 2.0),
                  CC = H0 / PetscPowReal(1.0 - 1.0 / n,qq);
  DMDALocalInfo info;
  DM            coordDA;
  Vec           coordinates;
  DMDACoor2d    **coords;
  PetscReal     **aH, r, s, tmp;
  PetscInt      j, k;

  PetscFunctionBeginUser;
  ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(user->da, &coordDA); CHKERRQ(ierr);
  ierr = DMGetCoordinates(user->da, &coordinates); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(coordDA, coordinates, &coords); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, H, &aH);CHKERRQ(ierr);
  for (k=info.ys; k<info.ys+info.ym; k++) {
      for (j=info.xs; j<info.xs+info.xm; j++) {
          r = radialcoord(coords[k][j]);
          if (r < L) {
              s = r / L;
              tmp = mm * s - (1.0/n) + PetscPowReal(1.0-s,mm) - PetscPowReal(s,mm);
              aH[k][j] = CC * PetscPowReal(tmp,qq);
          } else {
              aH[k][j] = 0.0;
          }
      }
  }
  ierr = DMDAVecRestoreArray(user->da, H, &aH);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode SetToSMB(Vec m, PetscBool chopneg, const AppCtx *user) {
  PetscErrorCode ierr;

  const PetscReal L  = user->exactL,
                  H0 = user->exactH0,
                  n  = user->n,
                  pp = 1.0 / n,
                  CC = user->Gamma * PetscPowReal(H0,2.0*n+2.0)
                          / PetscPowReal(2.0 * L * (1.0-1.0/n),n);
  DMDALocalInfo info;
  DM            coordDA;
  Vec           coordinates;
  DMDACoor2d    **coords;
  PetscReal     **am, r, s, tmp1, tmp2;
  PetscInt      j, k;

  PetscFunctionBeginUser;
  ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(user->da, &coordDA); CHKERRQ(ierr);
  ierr = DMGetCoordinates(user->da, &coordinates); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(coordDA, coordinates, &coords); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, m, &am);CHKERRQ(ierr);
  for (k=info.ys; k<info.ys+info.ym; k++) {
      for (j=info.xs; j<info.xs+info.xm; j++) {
          r = radialcoord(coords[k][j]);
          if (r > L - 0.01)  r = L - 0.01;
          s = r / L;
          tmp1 = PetscPowReal(s,pp) + PetscPowReal(1.0-s,pp) - 1.0;
          tmp2 = 2.0 * PetscPowReal(s,pp) + PetscPowReal(1.0-s,pp-1.0) * (1.0 - 2.0*s) - 1.0;
          am[k][j] = (CC / r) * PetscPowReal(tmp1,n-1.0) * tmp2;
          if (chopneg == PETSC_TRUE)  am[k][j] = PetscMax(am[k][j],0.0);
      }
  }
  ierr = DMDAVecRestoreArray(user->da, m, &am);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


// the single nontrivial formula for SIA
PetscReal getD(const PetscReal H, const PetscReal sx, const PetscReal sy, const AppCtx *user) {
    const PetscReal eps = user->eps,
                    n   = (1.0 - eps) * user->n + eps * 1.0;
    return user->Gamma * PetscPowReal(H,n+2.0) * PetscPowReal(sx*sx + sy*sy,(n-1.0)/2);
}


typedef struct {
    PetscReal x,y;
} Grad;


/* first of two nontrival operation with Q1 interpolants on an element */
/* evaluate the gradient of the surface elevation at any point (x,y) on element, using corner values of H and b */
PetscErrorCode gradsatpt(PetscInt j, PetscInt k,         // (j,k) is the element (by lower-left corner)
                         PetscReal locx, PetscReal locy, // = (x,y) coords in element
                         PetscReal **H, PetscReal **b,   // H[k][j] and b[k][j] are node values
                         const AppCtx *user, Grad *grads) {

  const PetscReal dx = user->dx, dy = dx,
                  x[4]  = {1.0 - locx / dx, locx / dx,       1.0 - locx / dx, locx / dx},
                  gx[4] = {- 1.0 / dx,      1.0 / dx,        - 1.0 / dx,      1.0 / dx},
                  y[4]  = {1.0 - locy / dy, 1.0 - locy / dy, locy / dy,       locy / dy},
                  gy[4] = {- 1.0 / dy,      - 1.0 / dy,      1.0 / dy,        1.0 / dy};
  PetscReal s[4];

  if ((grads == NULL) || (H == NULL) || (b == NULL)) {
      SETERRQ(PETSC_COMM_WORLD,1,"ERROR: illegal NULL ptr in gradsatpt() ...\n");
  }
  s[0] = H[k][j]     + b[k][j];
  s[1] = H[k][j+1]   + b[k][j+1];
  s[2] = H[k+1][j]   + b[k+1][j];
  s[3] = H[k+1][j+1] + b[k+1][j+1];
  grads->x = gx[0] * y[0] * s[0] + gx[1] * y[1] * s[1] + gx[2] * y[2] * s[2] + gx[3] * y[3] * s[3];
  grads->y =  x[0] *gy[0] * s[0] +  x[1] *gy[1] * s[1] +  x[2] *gy[2] * s[2] +  x[3] *gy[3] * s[3];
  PetscFunctionReturn(0);
}


Grad gradav(Grad g1, Grad g2) {
  Grad gav;
  gav.x = (g1.x + g2.x) / 2.0;
  gav.y = (g1.y + g2.y) / 2.0;
  return gav;
}


typedef struct {
    PetscReal x,y;
} Flux;


/* second of two nontrival operation with Q1 interpolants on an element */
/* evaluate the flux at any point (x,y) on element, using corner values of H and b */
PetscErrorCode fluxatpt(PetscInt j, PetscInt k,         // (j,k) is the element (by lower-left corner)
                        PetscReal locx, PetscReal locy, // = (x,y) coords in element
                        Grad gs,                        // compute with gradsatpt() first
                        PetscReal **H,                  // H[k][j] are node values
                        const AppCtx *user, Flux *q) {
  const PetscReal eps = user->eps,  dx = user->dx,  dy = dx,
                  x[4]  = {1.0 - locx / dx, locx / dx,       1.0 - locx / dx, locx / dx},
                  y[4]  = {1.0 - locy / dy, 1.0 - locy / dy, locy / dy,       locy / dy};
  PetscReal HH, DD;
  if ((q == NULL) || (H == NULL)) {
      SETERRQ(PETSC_COMM_WORLD,1,"ERROR: illegal NULL ptr in fluxatpt() ...\n");
  }
  HH = x[0] * y[0] * H[k][j] + x[1] * y[1] * H[k][j+1] + x[2] * y[2] * H[k+1][j] + x[3] * y[3] * H[k+1][j+1];
  DD = (1.0 - eps) * getD(HH,gs.x,gs.y,user) + eps * user->D0;
  q->x = - DD * gs.x;
  q->y = - DD * gs.y;
  PetscFunctionReturn(0);
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
  Grad            gsn, gss, gse, gsw;
  Flux            qn, qs, qe, qw;
  Vec             bloc, qloc;

  PetscFunctionBeginUser;
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
          if (user->mtrue == PETSC_FALSE) {  // default Mahaffy* method
              ierr = gradsatpt(j,k,     dx/2.0, 3.0*dy/4.0,      H,ab, user,&gsn); CHKERRQ(ierr);
              ierr =  fluxatpt(j,k,     dx/2.0, 3.0*dy/4.0, gsn, H,    user,&qn ); CHKERRQ(ierr);
              ierr = gradsatpt(j,k,     dx/2.0,     dy/4.0,      H,ab, user,&gss); CHKERRQ(ierr);
              ierr =  fluxatpt(j,k,     dx/2.0,     dy/4.0, gss, H,    user,&qs ); CHKERRQ(ierr);
              ierr = gradsatpt(j,k, 3.0*dx/4.0,     dy/2.0,      H,ab, user,&gse); CHKERRQ(ierr);
              ierr =  fluxatpt(j,k, 3.0*dx/4.0,     dy/2.0, gse, H,    user,&qe ); CHKERRQ(ierr);
              ierr = gradsatpt(j,k,     dx/4.0,     dy/2.0,      H,ab, user,&gsw); CHKERRQ(ierr);
              ierr =  fluxatpt(j,k,     dx/4.0,     dy/2.0, gsw, H,    user,&qw ); CHKERRQ(ierr);
          } else {  // true Mahaffy method
                    // this implementation is at least a factor of two inefficient
              Grad  gsn_nbr, gss_nbr, gse_nbr, gsw_nbr;
              // center-top point in element
              ierr =   gradsatpt(j,k,   dx/2.0, dy,   H,ab, user,&gsn); CHKERRQ(ierr);
              if (k < info->ys + info->ym - 1) {
                ierr = gradsatpt(j,k+1, dx/2.0, 0.0,  H,ab, user,&gsn_nbr); CHKERRQ(ierr);
                gsn  = gradav(gsn,gsn_nbr);
              }
              ierr =    fluxatpt(j,k,   dx/2.0, dy, gsn, H, user,&qn ); CHKERRQ(ierr);
              // center-bottom point in element
              ierr =   gradsatpt(j,k,   dx/2.0, 0.0,  H,ab, user,&gss); CHKERRQ(ierr);
              if (k > info->ys - 1) {
                ierr = gradsatpt(j,k-1, dx/2.0, dy,   H,ab, user,&gss_nbr); CHKERRQ(ierr);
                gss  = gradav(gss,gss_nbr);
              }
              ierr =    fluxatpt(j,k,   dx/2.0, 0.0,gss, H, user,&qs ); CHKERRQ(ierr);
              // center-right point in element
              ierr =   gradsatpt(j,k,   dx,  dy/2.0,  H,ab, user,&gse); CHKERRQ(ierr);
              if (j < info->xs + info->xm - 1) {
                ierr = gradsatpt(j+1,k, 0.0, dy/2.0,  H,ab, user,&gse_nbr); CHKERRQ(ierr);
                gse  = gradav(gse,gse_nbr);
              }
              ierr =    fluxatpt(j,k,   dx,  dy/2.0,gse, H, user,&qe ); CHKERRQ(ierr);
              // center-left point in element
              ierr =   gradsatpt(j,k,   0.0, dy/2.0,  H,ab, user,&gsw); CHKERRQ(ierr);
              if (j > info->xs - 1) {
                ierr = gradsatpt(j-1,k, dx,  dy/2.0,  H,ab, user,&gsw_nbr); CHKERRQ(ierr);
                gsw  = gradav(gsw,gsw_nbr);
              }
              ierr =    fluxatpt(j,k,   0.0, dy/2.0,gsw, H, user,&qw ); CHKERRQ(ierr);
          }
          aq[k][j] = (FluxQuad){qn.x,qs.x,qe.y,qw.y};
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

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"mah_","options to mahaffy","");CHKERRQ(ierr);
  ierr = PetscOptionsString(
      "-dump", "dump fields (H,Herror,m) into ascii files with this prefix",
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
      "-showdata", "use PETSc X viewers to show bed topography and climatic mass balance data",
      NULL,user->showdata,&user->showdata,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-true", "use true Mahaffy method, not default Mahaffy*",
      NULL,user->mtrue,&user->mtrue,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


// try a SNES solve; H returns modified
PetscErrorCode SNESAttempt(SNES snes, Vec H, PetscInt *its, SNESConvergedReason *reason) {
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  //ierr = VecScale(H,1.001); CHKERRQ(ierr);
  ierr = SNESSolve(snes, NULL, H); CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(snes,its);CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(snes,reason);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


// prints the |.|_inf = (max error) and |.|_1/(dim) = (av error) norms
PetscErrorCode ErrorReport(Vec H, DMDALocalInfo *info, AppCtx *user) {
  PetscErrorCode ierr;
  PetscScalar enorminf,enorm1;
  Vec         Hexact;
  if ((user->Ly != user->Lx) || (user->Lx != 900.0e3)) {
      SETERRQ(PETSC_COMM_WORLD,1,"exact thickness only defined for square 900km grid ...\n");
  }
  ierr = VecDuplicate(H,&Hexact); CHKERRQ(ierr);
  ierr = SetToExactThickness(Hexact,user); CHKERRQ(ierr);
  ierr = VecAXPY(Hexact,-1.0,H); CHKERRQ(ierr);    // Hexact < Hexact + (-1.0) H
  ierr = VecNorm(Hexact,NORM_INFINITY,&enorminf); CHKERRQ(ierr);
  ierr = VecNorm(Hexact,NORM_1,&enorm1); CHKERRQ(ierr);
  enorm1 /= info->mx * info->my;
  ierr = PetscPrintf(PETSC_COMM_WORLD,
             "     errors:  max = %8.2f,  av = %8.2f\n",enorminf,enorm1); CHKERRQ(ierr);
  ierr = VecDestroy(&Hexact);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


//  write a vector into a given filename
PetscErrorCode ViewToVTKASCII(Vec u, const char prefix[], const char name[]) {
    PetscErrorCode ierr;
    PetscViewer    viewer;
    int            strerr;
    char           filename[1024];

    strerr = sprintf(filename,"%s%s",prefix,name);
    if (strerr < 0) {
        SETERRQ1(PETSC_COMM_WORLD,6,"sprintf() returned %d < 0 ... stopping\n",strerr);
    }
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_VTK); CHKERRQ(ierr);
    ierr = VecView(u,viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


//  write various vectors to various files
PetscErrorCode DumpToFiles(Vec H, AppCtx *user) {
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    Vec            x, y;
    PetscInt       j, k;
    PetscReal      *ax, *ay;

    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF,info.mx,&x);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF,info.my,&y);CHKERRQ(ierr);
    ierr = VecGetArray(x, &ax);CHKERRQ(ierr);
    for (j=0; j<info.mx; j++) {
        ax[j] = -user->Lx + user->dx/2.0 + j * user->dx;
    }
    ierr = VecRestoreArray(x, &ax);CHKERRQ(ierr);
    ierr = VecGetArray(y, &ay);CHKERRQ(ierr);
    for (k=0; k<info.my; k++) {
        ay[k] = -user->Ly + user->dx/2.0 + k * user->dx;
    }
    ierr = VecRestoreArray(y, &ay);CHKERRQ(ierr);
    ierr = ViewToVTKASCII(x,user->figsprefix,"x.txt"); CHKERRQ(ierr);
    ierr = ViewToVTKASCII(y,user->figsprefix,"y.txt"); CHKERRQ(ierr);
    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = VecDestroy(&y); CHKERRQ(ierr);

    ierr = ViewToVTKASCII(H,user->figsprefix,"H.txt"); CHKERRQ(ierr);
    ierr = ViewToVTKASCII(user->b,user->figsprefix,"b.txt"); CHKERRQ(ierr);
    ierr = ViewToVTKASCII(user->m,user->figsprefix,"m.txt"); CHKERRQ(ierr);

    if (user->read == PETSC_FALSE) {
        Vec  Hexact;
        ierr = VecDuplicate(H,&Hexact); CHKERRQ(ierr);
        ierr = SetToExactThickness(Hexact,user); CHKERRQ(ierr);
        ierr = VecAXPY(Hexact,-1.0,H); CHKERRQ(ierr);    // do:  Hexact < Hexact + (-1.0) H
        ierr = VecScale(Hexact,-1.0); CHKERRQ(ierr);     // now: Hexact = H - Hexact
        ierr = ViewToVTKASCII(Hexact,user->figsprefix,"Herror.txt"); CHKERRQ(ierr);
        ierr = VecDestroy(&Hexact);CHKERRQ(ierr);
    } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,
             "NOTE: not writing Herror.txt file because grid and data were read from file\n"); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


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


PetscErrorCode ReadFromBinary(AppCtx *user) {
    PetscErrorCode ierr;
    PetscViewer viewer;
    Vec x, y;
    PetscReal *ax, *ay, fulllengthx, fulllengthy;

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"grn.dat",FILE_MODE_READ,&viewer); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&x); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)x,"x-axis from file"); CHKERRQ(ierr);
    ierr = VecLoad(x,viewer); CHKERRQ(ierr);
    // FIXME: check we own by: VecGetOwnershipRange(Vec x,PetscInt *low,PetscInt *high)
    ierr = VecGetSize(x,&user->Nx); CHKERRQ(ierr);
    if (user->Nx < 4) {  // 4 is somewhat arbitrary
        SETERRQ(PETSC_COMM_WORLD,1,"read Vec x has size too small\n");
    }
    ierr = VecGetArray(x, &ax);CHKERRQ(ierr);
    user->dx = ax[1] - ax[0];
    fulllengthx = ax[user->Nx-1] - ax[0] + user->dx;
    user->Lx = fulllengthx / 2.0;
    ierr = VecRestoreArray(x, &ax);CHKERRQ(ierr);
    ierr = VecDestroy(&x); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&y); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)y,"y-axis from file"); CHKERRQ(ierr);
    ierr = VecLoad(y,viewer); CHKERRQ(ierr);
    ierr = VecGetSize(y,&user->Ny); CHKERRQ(ierr);
    if (user->Ny < 4) {  // 4 is somewhat arbitrary
        SETERRQ(PETSC_COMM_WORLD,2,"read Vec y has size too small\n");
    }
    ierr = VecGetArray(y, &ay);CHKERRQ(ierr);
    fulllengthy = ay[user->Ny-1] - ay[0] + user->dx;
    user->Ly = fulllengthy / 2.0;
    ierr = VecRestoreArray(y, &ay);CHKERRQ(ierr);
    ierr = VecDestroy(&y); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&user->topgread); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)user->topgread,"bed topography from file"); CHKERRQ(ierr);
    ierr = VecLoad(user->topgread,viewer); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&user->cmbread); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)user->cmbread,"climatic mass balance from file"); CHKERRQ(ierr);
    ierr = VecLoad(user->cmbread,viewer); CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


PetscErrorCode ReshapeReadVecs(AppCtx *user) {
    PetscErrorCode ierr;
    PetscInt       j, k, s;
    PetscReal      *atopg, *acmb, **ab, **am;
    DMDALocalInfo  info;

    PetscFunctionBeginUser;
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = VecGetArray(user->topgread, &atopg);CHKERRQ(ierr);
    ierr = VecGetArray(user->cmbread, &acmb);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->b, &ab);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->m, &am);CHKERRQ(ierr);
    for (k=info.ys; k<info.ys+info.ym; k++) {
        for (j=info.xs; j<info.xs+info.xm; j++) {
            s = k * info.mx + j;  //  0 <= s <= info.my*info.mx - 1
            ab[k][j] = atopg[s];
            am[k][j] = acmb[s];
        }
    }
    ierr = VecRestoreArray(user->topgread, &atopg);CHKERRQ(ierr);
    ierr = VecRestoreArray(user->cmbread, &acmb);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->b, &ab);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->m, &am);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

