static const char help[] =
"Solves conservation-equation-for-layer problem in 1d:\n"
"    u_t + q_x = f\n"
"where the flux q combines an advecting layer\n"
"    q_0 = v(x) u\n"
"with a p-Laplacian\n"
"    q_1 = - k |u_x|^{p-2} u_x\n"
"so\n"
"    q = (1-lambda) q_0 + lambda q_1\n"
"and 0 <= lambda <= 1.\n"
"Domain is 0 < x < L, with periodic boundary conditions, subject to constraint\n"
"    u >= 0.\n"
"Uses SNESVI.  Several O(dx^2) finite difference methods to choose among, and\n"
"exact solution for v(x)=v0 case, and.\n\n";

//FIXME: add Jacobian

//FIXME: relate to event detection in \infty dimensions

//   ./layer -help |grep lay_

//   ./layer -snes_fd -draw_pause 0.5
//   ./layer -snes_mf -draw_pause 0.5

//   ./layer -snes_fd -lay_noshow
//   ./layer -snes_fd -lay_steps 100
//   ./layer -snes_fd -lay_dt 0.1
//   ./layer -snes_fd -lay_exactinit

//   ./layer -snes_fd -lay_scheme 1 -draw_pause 0.5
//   ./layer -snes_fd -lay_scheme 2 -draw_pause 0.5

//   ./layer -snes_fd -snes_type vinewtonssls
//   ./layer -snes_fd -snes_vi_monitor

//   ./layer -snes_fd -lay_dt 0.01 -lay_steps 1000 -da_refine 2

//   ./layer -snes_fd -lay_steps 1000 -da_refine 3 -lay_dt 0.0025 -snes_rtol 1.0e-10

// run at 10^5 CFL with 1.6 million DOFs
//   ./layer -lay_noshow -lay_steps 10 -da_refine 15 -lay_exactinit -lay_scheme 1

// for lev in 0 1 2 3 4 5 6; do ./layer -snes_fd -lay_exactinit -lay_noshow -lay_dt 0.01 -lay_steps 10 -lay_adscheme 1 -da_refine $lev | grep error; done

// ./layer -lay_steps 500 -lay_adscheme 2 -lay_lambda 0.9 -lay_k 100.0 -da_refine 3

#include <math.h>
#include <petscdmda.h>
#include <petscsnes.h>


typedef struct {
  DM        da;
  Vec       uold; // values of u[] at start of time step; u_{n-1}(x)
  PetscReal dt,   // fixed positive time step
            dx,   // fixed grid spacing
            L,    // domain is  0 <= x <= L
            v0,   // scale for velocity v(x)
            f0,   // scale for source term f(x)
            p,    // exponent for p-laplacian term
            k,    // scale for p-laplacian term
            lambda;    // 0.0 = advection, 1.0 = p-Laplacian
  PetscInt  adscheme,  // 0 = centered, 1 = third-order upwind-biased, 2 = box
            vchoice,   // 0 = (v(x) = v0), 1 = (v(x) quadratic)
            NN;        // number of time steps
  PetscBool exactinit; // initialize with exact solution
} AppCtx;


extern PetscErrorCode FormBounds(SNES,Vec,Vec);
extern PetscErrorCode SetToExactSolution(Vec,const AppCtx*);
extern PetscErrorCode FillVecs(Vec,Vec,Vec,const AppCtx*);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*,PetscScalar*,PetscScalar*,AppCtx*);
extern PetscErrorCode ProcessOptions(AppCtx*,PetscBool*,PetscBool*,char[]);
extern PetscErrorCode ViewToVTKASCII(Vec,const char[],const char[],const PetscInt);


int main(int argc,char **argv) {
  PetscErrorCode      ierr;
  SNES                snes;
  Vec                 u;
  AppCtx              user;
  PetscBool           noshow,genfigs;
  char                figsprefix[512];
  DMDALocalInfo       info;

  PetscInitialize(&argc,&argv,(char*)0,help);

  user.L  = 10.0;
  user.v0 = 10.0;
  user.f0 = 1.0;

  user.adscheme = 0;
  user.dt = 0.05;
  user.exactinit = PETSC_FALSE;
  user.k = 1.0;
  user.lambda = 0.0;
  user.p = 2.0;
  user.NN = 10;
  user.vchoice = 0;
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
  ierr = DMDASNESSetFunctionLocal(user.da,INSERT_VALUES,
            (PetscErrorCode (*)(DMDALocalInfo*,void*,void*,void*))FormFunctionLocal,
            &user);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr); CHKERRQ(ierr);

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
      if (reason < 0) {
          SETERRQ(PETSC_COMM_WORLD,3,"SNESVI solve diverged; stopping ...\n");
      }
      if (!noshow) {
          ierr = VecView(u, PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);
      }
      ierr = VecCopy(u, user.uold); CHKERRQ(ierr);
      t += user.dt;
    }
  }

  // evaluate numerical error relative to steady state
  if ((user.vchoice==0) && (user.lambda==0.0)) {
      Vec uexact;
      PetscScalar errnorm;
      ierr = DMCreateGlobalVector(user.da,&uexact);CHKERRQ(ierr);
      ierr = SetToExactSolution(uexact,&user);CHKERRQ(ierr);
      ierr = VecAXPY(u,-1.0,uexact); CHKERRQ(ierr);    // u <- u + (-1.0) uxact
      ierr = VecNorm(u,NORM_INFINITY,&errnorm); CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,
               "on dx=%.4e grid with N=%4d steps of dt=%.4e:  error |u-uexact|_inf = %g\n",
               user.dx,user.NN,user.dt,errnorm); CHKERRQ(ierr);
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
  if (user->vchoice > 0) {
    SETERRQ(PETSC_COMM_WORLD,4,
      "exact solution only available in v(x)=v0 case\n");
  }
  if (user->lambda > 0.0) {
    SETERRQ(PETSC_COMM_WORLD,5,
      "exact solution only available in pure-advection (lambda=0) case\n");
  }
  const PetscReal twopi = 2.0 * PETSC_PI,
                  x0    = (user->L/twopi) * asin(1.0/5.0),
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
        cosdiff = cos(twopi*x0/user->L) - cos(twopi*x/user->L);
        au[j]   = scale * ( (user->L/twopi) * cosdiff - (x-x0)/5.0 );
        au[j]   = PetscMax(0.0,au[j]);
      }
  }
  ierr = DMDAVecRestoreArray(user->da, u, &au);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


// constant velocity
PetscReal velocity(const PetscReal x, const AppCtx *user) {
  if (user->vchoice == 1) {
    const PetscReal Lsqr = user->L * user->L;
    return (4.0 / Lsqr) * user->v0 * x * (user->L - x);
  } else
    return user->v0;
}


// without constraint, with this f(x), \int_0^L u(t,x) dt --> - \infty
PetscReal fsource(const PetscReal x, const AppCtx *user) {
  const PetscReal CC = 2.0 * PETSC_PI / user-> L;
  return - (user->f0/5.0) + user->f0 * sin(CC * x);
}


// p-Laplacian
PetscReal plap(const PetscReal dudx, const PetscReal p) {
  return PetscPowReal(PetscAbsReal(dudx),p-2.0) * dudx;
}


/* for call-back: evaluate residual FF(x) on local process patch */
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info,PetscScalar *u,PetscScalar *FF,
                                 AppCtx *user) {
  PetscErrorCode  ierr;
  PetscInt        j;
  const PetscReal dt = user->dt, dx = user->dx, p = user->p,
                  nu = dt / dx, nu2 = nu / 2.0, nu4 = nu / 4.0;
  PetscReal       *uold;

  PetscFunctionBeginUser;
  ierr = DMDAVecGetArray(info->da, user->uold, &uold);CHKERRQ(ierr);
  for (j=info->xs; j<info->xs+info->xm; j++) {
      // non-flux part of residual
      const PetscReal x = dx/2 + dx * (PetscReal)j;
      FF[j] = u[j] - uold[j] - dt * fsource(x,user);
      // add p-laplacian flux part
      const PetscReal duright    = (u[j+1] - u[j]) / dx,
                      duleft     = (u[j] - u[j-1]) / dx,
                      duoldright = (uold[j+1] - uold[j]) / dx,
                      duoldleft  = (uold[j] - uold[j-1]) / dx;
      PetscReal plFF;
      plFF = (plap(duright,p) - plap(duleft,p)) - (plap(duoldright,p) - plap(duoldleft,p));
      FF[j] += user->lambda * (- user->k * nu2 * plFF);
      // add advection part
      PetscReal adFF;
      switch (user->adscheme) {
        case 0 : { // backward-Euler, centered scheme
          const PetscReal vleft  = velocity(x - dx/2,user),
                          vright = velocity(x + dx/2,user);
          adFF =   - nu2 * vleft * u[j-1]
                   + nu2 * (vright - vleft) * u[j]
                   + nu2 * vright * u[j+1];
          break; }
        case 1 : { // third-order upwind biased implicit thing; REQUIRES v(x)=v0>0
          if (user->vchoice > 0) {
            SETERRQ(PETSC_COMM_WORLD,4,"only v(x)=v0 is allowed with third-order scheme\n");
          }
          const PetscReal mu = user->v0 * nu / 6.0;
          adFF =        mu * u[j-2]
                  - 6.0*mu * u[j-1]
                  + 3.0*mu * u[j]
                  + 2.0*mu * u[j+1];
          break; }
        case 2 : { // trapezoid rule, centered scheme;
          const PetscReal vleft  = velocity(x - dx/2,user),
                          vright = velocity(x + dx/2,user);
          adFF =    nu4 * (vright * (u[j+1] + u[j]) - vleft * (u[j] + u[j-1]))
                  + nu4 * (vright * (uold[j+1] + uold[j]) - vleft * (uold[j] + uold[j-1]));
          break; }
        default :
          SETERRQ(PETSC_COMM_WORLD,2,"not allowed value of adscheme\n");
      }
      FF[j] += (1.0 - user->lambda) * adFF;
  }
  ierr = DMDAVecRestoreArray(info->da, user->uold, &uold);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode FillVecs(Vec vx, Vec vf, Vec vb, const AppCtx *user) {
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  const PetscReal dx = user->dx, pi = PETSC_PI;
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
      ab[j] = 1.0 * (cos(2.0*pi*x/10.0) + 1.0) - 0.5 * sin(2.0*pi*x/5.0);
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
      "-adscheme", "choose FD scheme for advection: 0 centered BEuler, 1 third-order, 2 centered trapezoid",
      NULL,user->adscheme,&user->adscheme,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-dt", "length of time step",
      NULL,user->dt,&user->dt,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-exactinit", "initialize with exact solution",
      NULL,user->exactinit,&user->exactinit,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString(
      "-genfigs", "generate one ascii file for each frame using this prefix",
      NULL,figsprefix,figsprefix,sizeof(figsprefix),genfigs); CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-k", "q_1 = - k |u_x|^{p-2} u_x  is p-laplacian flux; this sets k",
      NULL,user->k,&user->k,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-lambda", "q = (1-lambda) q_0 + lambda q_1 where q_0 is advective part and q_1 is p-laplacian",
      NULL,user->lambda,&user->lambda,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-noshow", "do not show solution with X",
      NULL,*noshow,noshow,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-p", "q_1 = - k |u_x|^{p-2} u_x  is p-laplacian flux; this sets p",
      NULL,user->p,&user->p,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt(
      "-steps", "number of time steps",
      NULL,user->NN,&user->NN,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt(
      "-vchoice","choose form of v(x): 0 constant, 1 quadratic",
      NULL,user->vchoice,&user->vchoice,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


//  we can write out the solution u and the source f into a given subdirectory (=prefix)
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

