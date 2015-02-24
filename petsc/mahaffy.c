static const char help[] =
"Solves steady ice sheet problem in 2d:\n"
"    div (q^x,q^y) = m,\n"
"    (q^x,q^y) = - Gamma H^{n+2} |grad s|^{n-1} grad s,\n"
"where  H(x,y)  is the ice thickness,  b(x,y)  is the bed elevation,\n"
"and  s(x,y) = H(x,y) + b(x,y).\n"
"Here  n > 1  and  Gamma = 2 A (rho g)^n / (n+2).\n"
"Domain is  -L < x < L,  -L < y < L,  with periodic boundary conditions.\n"
"Computed in flat bed case where analytical Jacobian is known.\n"
"Compares Mahaffy and Mahaffy* schemes.  Uses SNESVI.\n\n";

//   ./mahaffy -help |grep m_

#include <math.h>
#include <petscdmda.h>
#include <petscsnes.h>

typedef struct {
  DM        da;
  PetscReal dx,     // fixed grid spacing
            L,      // domain is [-L,L] x [-L,L]
            n,      // Glen exponent for SIA flux term
            g,      // acceleration of gravity
            rho,    // ice density
            secpera,// number of seconds in a year
            A,      // ice softness
            Gamma,  // coefficient for SIA flux term
            exactL, // radius of exact ice sheet
            exactH0;// center thickness of exact ice sheet
  PetscBool star;   // use Mahaffy* instead of Mahaffy
} AppCtx;


extern PetscErrorCode FormBounds(SNES,Vec,Vec);
extern PetscErrorCode SetToExactThickness(Vec,const AppCtx*);
extern PetscErrorCode SetToSMB(Vec,const AppCtx*);
extern PetscErrorCode ProcessOptions(AppCtx*);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*,PetscScalar*,PetscScalar*,AppCtx*);
//see layer.c for this:
//  extern PetscErrorCode FormJacobianLocal(DMDALocalInfo*,PetscScalar*,Mat,Mat,AppCtx*);


int main(int argc,char **argv) {
  PetscErrorCode      ierr;
  SNES                snes;
  Vec                 H;
  AppCtx              user;
  DMDALocalInfo       info;

  PetscInitialize(&argc,&argv,(char*)0,help);

  user.L      = 900.0e3;    // m
  user.n      = 3.0;
  user.g      = 9.81;       // m/s^2
  user.rho    = 910.0;      // kg/m^3
  user.secpera= 31556926.0;
  user.A      = 1.0E-16/user.secpera; // = 3.17e-24  1/(Pa^3 s); EISMINT I value
  user.Gamma  = 2.0 * PetscPowReal(user.rho*user.g,user.n) * user.A / (user.n+2.0);

  user.exactL = 750.0e3;    // m
  user.exactH0= 3600.0;     // m

  ierr = ProcessOptions(&user); CHKERRQ(ierr);

  ierr = DMDACreate2d(PETSC_COMM_WORLD,
                      DM_BOUNDARY_PERIODIC,DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX,
                      -18,-18,PETSC_DECIDE,PETSC_DECIDE,  // initial 18x18 grid has dx=100km
                      1, 1,                               // dof=1,  stencilwidth=1
                      NULL,NULL,&user.da);
  ierr = DMSetApplicationContext(user.da, &user);CHKERRQ(ierr);

  ierr = DMDAGetLocalInfo(user.da,&info); CHKERRQ(ierr);
  user.dx = 2 * user.L / (PetscReal)(info.mx);
  if (info.mx == info.my)  {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "solving on %d x %d grid spacing dx=%g ...\n",
                         info.mx,info.my,user.dx); CHKERRQ(ierr);
  } else {
      SETERRQ(PETSC_COMM_WORLD,1,"mx != my ERROR: grid must have equal spacing in x,y ...\n");
  }

  // cell-centered grid
  ierr = DMDASetUniformCoordinates(user.da,
             -user.L+user.dx/2,user.L-user.dx/2,
             -user.L+user.dx/2,user.L-user.dx/2,
             0.0,1.0); CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(user.da,&H);CHKERRQ(ierr);
  ierr = VecSetOptionsPrefix(H,"H_"); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)H,"thickness solution H"); CHKERRQ(ierr);

  //ierr = SetToExactThickness(H,&user);CHKERRQ(ierr);
  ierr = VecSet(H,0.0); CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetDM(snes,user.da);CHKERRQ(ierr);
  ierr = SNESSetType(snes,SNESVINEWTONRSLS);CHKERRQ(ierr);
  ierr = SNESVISetComputeVariableBounds(snes,&FormBounds);CHKERRQ(ierr);
  ierr = DMDASNESSetFunctionLocal(user.da,INSERT_VALUES,
              (DMDASNESFunction)FormFunctionLocal,&user); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  PetscInt            its;
  SNESConvergedReason reason;
  // ierr = VecScale(H,0.99); CHKERRQ(ierr); // move so that line search goes some distance
  ierr = SNESSolve(snes, NULL, H); CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%3d Newton iterations (%s)\n",
                     its,SNESConvergedReasons[reason]);CHKERRQ(ierr);
  if (reason < 0) {
      SETERRQ(PETSC_COMM_WORLD,2,"SNESVI solve diverged; stopping ...\n");
  }

  PetscScalar errnorm;
  Vec         Hexact;
  ierr = VecDuplicate(H,&Hexact); CHKERRQ(ierr);
  ierr = SetToExactThickness(Hexact,&user); CHKERRQ(ierr);
  ierr = VecAXPY(H,-1.0,Hexact); CHKERRQ(ierr);    // H < - H + (-1.0) Hexact
  ierr = VecNorm(H,NORM_INFINITY,&errnorm); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"error |u-uexact|_inf = %g\n",errnorm); CHKERRQ(ierr);
  ierr = VecDestroy(&Hexact);CHKERRQ(ierr);

  ierr = VecDestroy(&H);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = DMDestroy(&user.da);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}


//  for call-back: tell SNESVI (variational inequality) that we want  0.0 <= H < +infinity
PetscErrorCode FormBounds(SNES snes, Vec Xl, Vec Xu) {
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = VecSet(Xl,0.0); CHKERRQ(ierr);
  ierr = VecSet(Xu,PETSC_INFINITY); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode SetToExactThickness(Vec H, const AppCtx *user) {
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  const PetscReal L  = user->exactL,
                  H0 = user->exactH0,
                  n  = user->n,
                  mm = 1.0 + 1.0 / n,
                  qq = n / (2.0 * n + 2.0);
  DMDALocalInfo info;
  PetscReal     *aH;
  PetscInt      j, k;
  ierr = DMDAGetLocalInfo(user->da,&info); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, H, &aH);CHKERRQ(ierr);
  for (j=info.xs; j<info.xs+info.xm; j++) {
      for (k=info.ys; j<info.ys+info.ym; j++) {
          r = FIXME;
          if (r < L) {
              PetscReal s, lamhat, Hs;
              if (r < 0.01)
                  r = 0.01;  // avoid r=0 singularity
              s = r / L;
              lamhat = mm * s - (1.0/n) + pow(1-s,mm) - pow(s,mm);
              Hs = (H0 / pow(1-1/n,qq)) * pow(lamhat,qq);
              aH[j][k] = Hs;
          } else {
              aH[j][k] = 0.0;
          }
      }
  }
  ierr = DMDAVecRestoreArray(user->da, H, &aH);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscReal SetToSMB(Vec m, const AppCtx *user) {
  FIXME
  return FIXME;
}


// flux formula for SIA
PetscErrorCode q(const PetscReal H, const PetscReal sx, const PetscReal sy, const AppCtx *user,
                 PetscReal *qx, PetscReal *qy) {
  PetscReal n = user->n, D;
  D = user->Gamma * PetscPowReal(H,n+2.0) * PetscPowReal(sx*sx + sy*sy,(n-1.0)/2);
  *qx = - D * sx;
  *qy = - D * sy;
  PetscFunctionReturn(0);
}


/* for call-back: evaluate residual FF(x) on local process patch */
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info,PetscScalar *u,PetscScalar *FF,
                                 AppCtx *user) {
  PetscErrorCode  ierr;
  const PetscReal dt = user->dt, dx = user->dx,
                  lam = user->lambda, v0 = user->v0,
                  nu = dt / dx;
  PetscInt        j;
  PetscReal       *uold;
  Vec             uoldloc;

  PetscFunctionBeginUser;
  ierr = DMCreateLocalVector(info->da,&uoldloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(info->da,user->uold,INSERT_VALUES,uoldloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(info->da,user->uold,INSERT_VALUES,uoldloc);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(info->da, uoldloc, &uold);CHKERRQ(ierr);
  for (j=info->xs; j<info->xs+info->xm; j++) {
      // non-flux part of residual
      const PetscReal x = dx/2 + dx * (PetscReal)j;
      FF[j] = u[j] - uold[j] - dt * fsource(user->midtime,x,user);
      // add p-laplacian (SIA) flux part q^0
      const PetscReal
          dbright    = (bedelevation(x+dx,user) - bedelevation(x,user)) / dx,
          dbleft     = (bedelevation(x,user) - bedelevation(x-dx,user)) / dx;
      const PetscReal
          uright     = (u[j] + u[j+1]) / 2.0,
          uleft      = (u[j-1] + u[j]) / 2.0,
          duright    = (u[j+1] - u[j]) / dx,
          duleft     = (u[j] - u[j-1]) / dx;
      const PetscReal
          uoldright  = (uold[j] + uold[j+1]) / 2.0,
          uoldleft   = (uold[j-1] + uold[j]) / 2.0,
          duoldright = (uold[j+1] - uold[j]) / dx,
          duoldleft  = (uold[j] - uold[j-1]) / dx;
      FF[j] += lam * (nu / 2.0) *
               (  S(uright,   duright   +dbright,user) - S(uleft,   duleft   +dbleft,user)
                + S(uoldright,duoldright+dbright,user) - S(uoldleft,duoldleft+dbleft,user) );
      // add advection part q^1
      PetscReal adFF, c[4];
      ierr = setadconstants(user->adscheme,c); CHKERRQ(ierr);
      adFF = v0 * (nu / 2.0) * ( c[0] * u[j-2] + c[1] * u[j-1] + c[2] * u[j] + c[3] * u[j+1] );
      if (user->adscheme == 2)
          adFF += v0 * (nu / 2.0) * (1.0 / 2.0) * ( uold[j+1] - uold[j-1] );
      FF[j] += (1.0 - lam) * adFF;
  }
  ierr = DMDAVecRestoreArray(info->da, uoldloc, &uold);CHKERRQ(ierr);

  ierr = VecDestroy(&uoldloc); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/* for call-back: evaluate residual FF(x) on local process patch */
PetscErrorCode ProcessOptions(AppCtx *user, PetscBool *noshow, PetscBool *genfigs, char figsprefix[], PetscBool *genmassfile, char massfilename[]) {
  PetscErrorCode ierr;

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"lay_","options to layer","");CHKERRQ(ierr);
  ierr = PetscOptionsInt(
      "-adscheme", "choose FD scheme for advection: 0 centered BEuler, 1 third-order BEuler, 2 centered trapezoid",
      NULL,user->adscheme,&user->adscheme,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-dt", "length of time step",
      NULL,user->dt,&user->dt,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-exactinit", "initialize with exact solution",
      NULL,user->exactinit,&user->exactinit,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-timedependentsource", "use an f(t,x) formula which is actually time-dependent",
      NULL,user->fusest,&user->fusest,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-gamma", "q_1 = - gamma u^{n+2} |(u+b)_x|^{n-1} (u+b)_x  is SIA flux; this sets gamma",
      NULL,user->gamma,&user->gamma,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString(
      "-genfigs", "if set, generate one ascii file for each frame using this prefix",
      NULL,figsprefix,figsprefix,512,genfigs); CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-glenn", "q_1 = - gamma u^{n+2} |(u+b)_x|^{n-1} (u+b)_x  is SIA flux; this sets n",
      NULL,user->n,&user->n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-jac", "use analytical jacobian; default is to do finite difference (-snes_fd) or matrix-free (-snes_mf)",
      NULL,user->usejacobian,&user->usejacobian,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-lambda", "q = lambda q^0 + (1-lambda) q^1 where q^0 is SIA part and q^1 is advective part",
      NULL,user->lambda,&user->lambda,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString(
      "-massfile", "if set, write mass time-series [t_n M_n R_n C_n balance] to this file",
      NULL,massfilename,massfilename,512,genmassfile);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-noshow", "do _not_ show solution with X window viewers",
      NULL,*noshow,noshow,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt(
      "-steps", "number of time steps",
      NULL,user->NN,&user->NN,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}



