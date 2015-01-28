static const char help[] =
"Solves conservation-equation-for-layer problem in 1d:\n"
"    u_t + q_x = f,\n"
"where the flux\n"
"    q = lambda q^0 + (1-lambda) q^1,\n"
"with  0 <= lambda <= 1, combines an SIA-type flux\n"
"    q^0 = - gamma u^{n+2} |(u+b)_x|^{n-1} (u+b)_x\n"
"with an advecting layer\n"
"    q^1 = v0 u.\n"
"Here n >= 1 and b(x) is a smooth function given in bedelevation() below.\n"
"Domain is 0 < x < L, with periodic boundary conditions, subject to constraint\n"
"    u >= 0.\n"
"Uses SNESVI.  Several O(dx^2) finite difference methods to choose among.\n"
"Exact solution for lambda=0 case.  Either analytical Jacobian or\n"
"finite-difference evaluation of Jacobian.\n\n";

//FIXME: relate to event detection in \infty dimensions

//   ./layer -help |grep lay_

//   ./convtest.sh

//   ./layer -lay_jac

//   ./layer -snes_mf
//   ./layer -snes_type vinewtonssls
//   ./layer -lay_noshow -snes_vi_monitor -snes_rtol 1.0e-10

//   ./layer -lay_dt 0.1 -lay_exactinit -da_refine 3 -lay_adscheme 1

// ./layer -lay_steps 500 -draw_pause 0.0 -lay_adscheme 1 -da_refine 3 -lay_lambda 1.0 -lay_dt 0.02 -lay_gamma 0.1 -lay_glenn 3.0

// run steady state at 10^5 CFL with 1.6 million DOFs (fd Jacobian)
//   ./layer -lay_noshow -lay_steps 5 -da_refine 15 -lay_exactinit -lay_adscheme 1

// run steady state at 6 x 10^5 CFL with 13 million DOFs (analytical Jacobian)
//   ./layer -lay_noshow -lay_steps 5 -da_refine 18 -lay_exactinit -lay_adscheme 1 -lay_jac


#include <math.h>
#include <petscdmda.h>
#include <petscsnes.h>


typedef struct {
  DM        da;
  Vec       uold; // values of u[] at start of time step; u_{n-1}(x)
  PetscReal dt,   // fixed positive time step
            dx,   // fixed grid spacing
            L,    // domain is  0 <= x <= L
            v0,   // velocity (constant)
            f0,   // scale for source term f(x)
            n,    // Glen exponent for SIA flux term
            gamma,// coefficient for SIA flux term
            lambda;    // 0.0 = advection, 1.0 = SIA
  PetscInt  adscheme,  // 0 = centered, 1 = third-order upwind-biased, 2 = box
            NN;        // number of time steps
  PetscBool exactinit, // initialize with exact solution
            usejacobian;// use analytical jacobian
} AppCtx;


extern PetscErrorCode FormBounds(SNES,Vec,Vec);
extern PetscErrorCode SetToExactSolution(Vec,const AppCtx*);
extern PetscErrorCode FillVecs(Vec,Vec,Vec,const AppCtx*);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*,PetscScalar*,PetscScalar*,AppCtx*);
extern PetscErrorCode FormJacobianLocal(DMDALocalInfo*,PetscScalar*,Mat,Mat,AppCtx*);
extern PetscErrorCode ProcessOptions(AppCtx*,PetscBool*,PetscBool*,char[]);
extern PetscErrorCode ViewToVTKASCII(Vec,const char[],const char[],const PetscInt);


int main(int argc,char **argv) {
  PetscErrorCode      ierr;
  SNES                snes;
  Vec                 u;
  AppCtx              user;
  PetscBool           noshow,genfigs;
  char                figsprefix[512] = "";
  DMDALocalInfo       info;

  PetscInitialize(&argc,&argv,(char*)0,help);

  user.L  = 10.0;
  user.v0 = 10.0;
  user.f0 = 1.0;

  user.adscheme = 0;
  user.dt = 0.05;
  user.exactinit = PETSC_FALSE;
  user.usejacobian = PETSC_FALSE;
  user.n      = 3.0;
  user.gamma  = 1.0;
  user.lambda = 0.0;
  user.NN = 10;
  noshow = PETSC_FALSE;
  ierr = ProcessOptions(&user,&noshow,&genfigs,figsprefix); CHKERRQ(ierr);

  ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC,
                      -50,         // override with -da_grid_x or -da_refine
                      1, 2, NULL,  // dof = 1 and stencil width = 2
                      &user.da); CHKERRQ(ierr);
  ierr = DMSetApplicationContext(user.da, &user);CHKERRQ(ierr);

  ierr = DMDAGetLocalInfo(user.da,&info); CHKERRQ(ierr);
  user.dx = user.L / (PetscReal)(info.mx);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "doing %d steps of dt=%f on grid of %d points with spacing dx=%g\n"
                     "[CFL time step is %g]\n",
                     user.NN,user.dt,info.mx,user.dx,user.dx/user.v0); CHKERRQ(ierr);
  // cell-centered grid
  ierr = DMDASetUniformCoordinates(user.da,user.dx/2,user.L-user.dx/2,
                                   -1.0,-1.0,-1.0,-1.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(user.da,&u);CHKERRQ(ierr);
  ierr = VecSetOptionsPrefix(u,"u_"); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)u,"solution u"); CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(user.da,&user.uold);CHKERRQ(ierr);
  ierr = VecSetOptionsPrefix(user.uold,"uold_"); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)user.uold,"OLD solution uold"); CHKERRQ(ierr);

  if (user.exactinit) {
    ierr = SetToExactSolution(user.uold,&user);CHKERRQ(ierr);
  } else {
    ierr = VecSet(user.uold,0.0); CHKERRQ(ierr);
  }

  if (genfigs) {
    Vec x,f,b;
    ierr = VecDuplicate(u,&x); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)x,"source f"); CHKERRQ(ierr);
    ierr = VecDuplicate(u,&f); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)f,"coordinate x"); CHKERRQ(ierr);
    ierr = VecDuplicate(u,&b); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)b,"bed elevation b"); CHKERRQ(ierr);
    ierr = FillVecs(x,f,b,&user); CHKERRQ(ierr);
    ierr = ViewToVTKASCII(x,figsprefix,"x",-1); CHKERRQ(ierr);
    ierr = ViewToVTKASCII(f,figsprefix,"f",-1); CHKERRQ(ierr);
    ierr = ViewToVTKASCII(b,figsprefix,"b",-1); CHKERRQ(ierr);
    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = VecDestroy(&f); CHKERRQ(ierr);
    ierr = VecDestroy(&b); CHKERRQ(ierr);
  }

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetDM(snes,user.da);CHKERRQ(ierr);
  ierr = SNESSetType(snes,SNESVINEWTONRSLS);CHKERRQ(ierr);
  ierr = SNESVISetComputeVariableBounds(snes,&FormBounds);CHKERRQ(ierr);
  ierr = DMDASNESSetFunctionLocal(user.da,INSERT_VALUES,(DMDASNESFunction)FormFunctionLocal,&user); CHKERRQ(ierr);
  if (user.usejacobian) {
    ierr = DMDASNESSetJacobianLocal(user.da,(DMDASNESJacobian)FormJacobianLocal,&user); CHKERRQ(ierr);
  }
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  /* time-stepping loop */
  {
    PetscReal  t;
    PetscInt   n, its;
    SNESConvergedReason reason;
    t = 0.0;
    ierr = VecCopy(user.uold,u);CHKERRQ(ierr);
    if (!noshow) {
      ierr = VecView(u, PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);
    }
    for (n = 0; n < user.NN; ++n) {
      ierr = PetscPrintf(PETSC_COMM_WORLD, "  time[%3d]=%6g: ", n+1, t+user.dt); CHKERRQ(ierr);
      ierr = VecScale(u,0.99); CHKERRQ(ierr); // move u so that line search goes some distance
      ierr = SNESSolve(snes, NULL, u); CHKERRQ(ierr);
      ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
      ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"%3d Newton iterations (%s)\n",
                         its,SNESConvergedReasons[reason]);CHKERRQ(ierr);
      if (genfigs) {
          ierr = ViewToVTKASCII(u,figsprefix,"u",n+1); CHKERRQ(ierr);
      }
      if (!noshow) {
          ierr = VecView(u, PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);
      }
      if (reason < 0) {
          SETERRQ(PETSC_COMM_WORLD,3,"SNESVI solve diverged; stopping ...\n");
      }
      ierr = VecCopy(u, user.uold); CHKERRQ(ierr);
      t += user.dt;
    }
  }

  // evaluate numerical error relative to steady state
  if (user.lambda==0.0) {
      Vec uexact;
      PetscScalar errnorm;
      ierr = DMCreateGlobalVector(user.da,&uexact);CHKERRQ(ierr);
      ierr = SetToExactSolution(uexact,&user);CHKERRQ(ierr);
      ierr = VecAXPY(u,-1.0,uexact); CHKERRQ(ierr);    // u <- u + (-1.0) uxact
      ierr = VecNorm(u,NORM_INFINITY,&errnorm); CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,
               "on dx=%.4e grid (%6d pts) with N=%4d steps of dt=%.4e:  error |u-uexact|_inf = %g\n",
               user.dx,info.mx,user.NN,user.dt,errnorm); CHKERRQ(ierr);
  }

  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&user.uold);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = DMDestroy(&user.da);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}


//  for call-back: tell SNESVI (variational inequality) that we want  0.0 <= u < +infinity
PetscErrorCode FormBounds(SNES snes, Vec Xl, Vec Xu) {
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = VecSet(Xl,0.0); CHKERRQ(ierr);
  ierr = VecSet(Xu,PETSC_INFINITY); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode SetToExactSolution(Vec u, const AppCtx *user) {
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  if (user->lambda > 0.0) {
    SETERRQ(PETSC_COMM_WORLD,5,
      "exact solution only available in pure-advection (lambda=0) case\n");
  }
  const PetscReal twopi = 2.0 * PETSC_PI,
                  L     = user->L,
                  shift = L / 15.0,
                  fdown = - 1.0 / 5.0;
  const PetscReal
                  x0    = (L/twopi) * asin(-fdown) + shift,
                  scale = user->f0/user->v0,
                  dx    = user->dx;
  DMDALocalInfo info;
  PetscReal     *au, x, cosdiff;
  PetscInt      j;
  ierr = DMDAGetLocalInfo(user->da,&info); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, u, &au);CHKERRQ(ierr);
  for (j=info.xs; j<info.xs+info.xm; j++) {
      x = dx/2 + dx * (PetscReal)j;
      if (x <= x0) {
        au[j] = 0.0;
      } else {
        cosdiff = cos((twopi/L)*(x0-shift)) - cos((twopi/L)*(x-shift));
        au[j]   = scale * ( (L/twopi) * cosdiff + fdown * (x-x0) );
        au[j]   = PetscMax(0.0,au[j]);
      }
  }
  ierr = DMDAVecRestoreArray(user->da, u, &au);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


// without constraint, with this f(x), \int_0^L u(t,x) dt --> - \infty
PetscReal fsource(const PetscReal x, const AppCtx *user) {
  const PetscReal twopi = 2.0 * PETSC_PI,
                  L     = user->L,
                  shift = L / 15.0,
                  fdown = - 1.0 / 5.0;
  return user->f0 * ( fdown + sin((twopi/L) * (x - shift)) );
}


// bed elevation for SIA-type flux
PetscReal bedelevation(const PetscReal x, const AppCtx *user) {
  const PetscReal pi = PETSC_PI, L = user->L, L2 = L / 2.0,
                  xx = x - user->L/15.0;
  return 1.0 * (cos(2.0*pi*xx/L) + 1.0) - 0.5 * sin(2.0*pi*xx/L2);
}


// flux formula for SIA, and its partial derivatives
PetscReal S(const PetscReal u, const PetscReal dhdx, const AppCtx *user) {
  return -user->gamma * PetscPowReal(u,user->n+2.0) * PetscPowReal(PetscAbsReal(dhdx),user->n-1.0) * dhdx;
}

PetscReal S1(const PetscReal u, const PetscReal dhdx, const AppCtx *user) {
  return -user->gamma * (user->n+2.0) * PetscPowReal(u,user->n+1.0) * PetscPowReal(PetscAbsReal(dhdx),user->n-1.0) * dhdx;
}

PetscReal S2(const PetscReal u, const PetscReal dhdx, const AppCtx *user) {
  return -user->gamma * PetscPowReal(u,user->n+2.0) * (user->n) * PetscPowReal(PetscAbsReal(dhdx),user->n-1.0);
}


PetscErrorCode setadconstants(PetscInt scheme, PetscReal *c) {
  switch (scheme) {
    case 0 : { // backward-Euler, centered scheme
      c[0] = 0.0;
      c[1] = -1.0/2.0;
      c[2] = 0.0;
      c[3] = 1.0/2.0;
      break; }
    case 1 : { // third-order upwind biased backward-Euler
      c[0] = 1.0/6.0;
      c[1] = -6.0/6.0;
      c[2] = 3.0/6.0;
      c[3] = 2.0/6.0;
      break; }
    case 2 : { // trapezoid rule, centered scheme;
      c[0] = 0.0;
      c[1] = -1.0/4.0;
      c[2] = 0.0;
      c[3] = 1.0/4.0;
      break; }
    default :
      SETERRQ(PETSC_COMM_WORLD,1,"not allowed value of scheme\n");
  }
  return 0;
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
      FF[j] = u[j] - uold[j] - dt * fsource(x,user);
      // add p-laplacian (SIA) flux part q^0
      const PetscReal
          dbright    = (bedelevation(x+dx,user) - bedelevation(x,user)) / dx,
          dbleft     = (bedelevation(x,user) - bedelevation(x-dx,user)) / dx;
      const PetscReal
          uright     = (u[j] + u[j+1]) / 2.0,
          uleft      = (u[j-1] + u[j]) / 2.0,
          duright    = (u[j+1] - u[j]) / dx,
          duleft     = (u[j] - u[j-1]) / dx,
          dhright    = duright + dbright,
          dhleft     = duleft + dbleft,
          dq         = S(uright,dhright,user) - S(uleft,dhleft,user);
      const PetscReal
          uoldright  = (uold[j] + uold[j+1]) / 2.0,
          uoldleft   = (uold[j-1] + uold[j]) / 2.0,
          duoldright = (uold[j+1] - uold[j]) / dx,
          duoldleft  = (uold[j] - uold[j-1]) / dx,
          dqold      = S(uoldright,duoldright+dbright,user) - S(uoldleft,duoldleft+dbleft,user);
      FF[j] += lam * (nu / 2.0) * (dq + dqold);
      // add advection part q^1
      PetscReal adFF = 0.0, c[4];
      ierr = setadconstants(user->adscheme,c); CHKERRQ(ierr);
      if (user->adscheme == 2)
          adFF += v0 * (nu / 4.0) * ( uold[j+1] - uold[j-1] );
      adFF += v0 * nu * ( c[0] * u[j-2] + c[1] * u[j-1] + c[2] * u[j] + c[3] * u[j+1] );
      FF[j] += (1.0 - lam) * adFF;
  }
  ierr = DMDAVecRestoreArray(info->da, uoldloc, &uold);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/* for call-back: evaluate Jacobian, for rows corresponding to local process patch */
PetscErrorCode FormJacobianLocal(DMDALocalInfo *info, PetscScalar *u, Mat jacpre, Mat jac,
                                 AppCtx *user) {
  PetscErrorCode  ierr;
  const PetscReal dt = user->dt,  dx = user->dx,
                  lam = user->lambda, v0 = user->v0,
                  nu = dt / dx;
  PetscInt        j, k;
  PetscReal       v[4];
  MatStencil      row,col[4];

  PetscFunctionBeginUser;
  for (j=info->xs; j<info->xs+info->xm; j++) {
      row.i = j;

      const PetscReal x = dx/2 + dx * (PetscReal)j;
      const PetscReal
          dbright    = (bedelevation(x+dx,user) - bedelevation(x,user)) / dx,
          dbleft     = (bedelevation(x,user) - bedelevation(x-dx,user)) / dx;
      const PetscReal
          uright     = (u[j] + u[j+1]) / 2.0,
          uleft      = (u[j-1] + u[j]) / 2.0,
          duright    = (u[j+1] - u[j]) / dx,
          duleft     = (u[j] - u[j-1]) / dx,
          dhright    = duright + dbright,
          dhleft     = duleft + dbleft;
      const PetscReal
          S1left     = S1(uleft, dhleft, user),
          S2left     = S2(uleft, dhleft, user),
          S1right    = S1(uright, dhright, user),
          S2right    = S2(uright, dhright, user);

      col[0].i = j-2;
      v[0] = 0.0;
      col[1].i = j-1;
      v[1] = lam * (nu / 2.0) * ( - S1left * (1.0/2.0) - S2left * (-1.0/dx) );
      col[2].i = j;
      v[2] = 1.0
             + lam * (nu / 2.0) * ( S1right * (1.0/2.0) + S2right * (-1.0/dx)
                                   - S1left * (1.0/2.0) - S2left * (1.0/dx)   );
      col[3].i = j+1;
      v[3] = lam * (nu / 2.0) * ( S1right * (1.0/2.0) + S2right * (1.0/dx) );

      PetscReal c[4];
      ierr = setadconstants(user->adscheme,c); CHKERRQ(ierr);
      for (k=0; k<4; k++)
          v[k] += (1.0 - lam) * v0 * nu * c[k];

      ierr = MatSetValuesStencil(jac,1,&row,4,col,v,INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatSetOption(jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode FillVecs(Vec vx, Vec vf, Vec vb, const AppCtx *user) {
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  const PetscReal dx = user->dx;
  DMDALocalInfo info;
  PetscReal     x, *ax, *af, *ab;
  PetscInt      j;
  ierr = DMDAGetLocalInfo(user->da,&info); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, vx, &ax);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, vf, &af);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, vb, &ab);CHKERRQ(ierr);
  for (j=info.xs; j<info.xs+info.xm; j++) {
      x     = dx/2 + dx * (PetscReal)j;
      ax[j] = x;
      af[j] = fsource(x,user);
      ab[j] = bedelevation(x,user);
  }
  ierr = DMDAVecRestoreArray(user->da, vx, &ax);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, vf, &af);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, vb, &ab);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/* for call-back: evaluate residual FF(x) on local process patch */
PetscErrorCode ProcessOptions(AppCtx *user, PetscBool *noshow, PetscBool *genfigs, char figsprefix[]) {
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
  ierr = PetscOptionsReal(
      "-gamma", "q_1 = - gamma u^{n+2} |(u+b)_x|^{n-1} (u+b)_x  is SIA flux; this sets gamma",
      NULL,user->gamma,&user->gamma,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString(
      "-genfigs", "generate one ascii file for each frame using this prefix",
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
  ierr = PetscOptionsBool(
      "-noshow", "do _not_ show solution with X window viewers",
      NULL,*noshow,noshow,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt(
      "-steps", "number of time steps",
      NULL,user->NN,&user->NN,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


//  write a vector into a given subdirectory (=prefix), with filename
//      prefix/name.txt
//  or
//      prefix/name-nn.txt    if nn >= 0
PetscErrorCode ViewToVTKASCII(Vec u, const char prefix[], const char name[],
                              const PetscInt nn) {
    PetscErrorCode ierr;
    PetscViewer viewer;
    char filename[1024];
    int  strerr;
    if (nn >= 0)
        strerr = sprintf(filename,"%s%s-%d.txt",prefix,name,nn);
    else
        strerr = sprintf(filename,"%s%s.txt",prefix,name);
    if (strerr < 0) {
        SETERRQ1(PETSC_COMM_WORLD,6,"sprintf() returned %d < 0 ... stopping\n",strerr);
    }
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_VTK); CHKERRQ(ierr);
    ierr = VecView(u,viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

