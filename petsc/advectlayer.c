static const char help[] = "Solves advecting-layer problem in 1d:\n"
"    u_t + div q = f\n"
"where\n"
"    q = v(x) u\n"
"on 0 < x < L, with periodic boundary conditions, subject to constraint\n"
"    u >= 0.\n"
"Implicit (backward Euler) centered finite difference method and\n"
"SNESVI.\n\n";

//FIXME: add Jacobian

//FIXME: relate to event detection in \infty dimensions

//   ./advectlayer -help |grep al_

//   ./advectlayer -snes_fd -draw_pause 0.5
//   ./advectlayer -snes_mf -draw_pause 0.5

//   ./advectlayer -snes_fd -al_noshow
//   ./advectlayer -snes_fd -al_steps 100
//   ./advectlayer -snes_fd -al_dt 0.1
//   ./advectlayer -snes_fd -al_exactinit

//   ./advectlayer -snes_fd -al_thirdorder -draw_pause 0.5
//   ./advectlayer -snes_fd -al_box  -draw_pause 0.5

//   ./advectlayer -snes_fd -snes_type vinewtonssls
//   ./advectlayer -snes_fd -snes_vi_monitor

//   ./advectlayer -snes_fd -al_dt 0.01 -al_steps 1000 -da_refine 2

//   ./advectlayer -snes_fd -al_steps 1000 -da_refine 3 -al_dt 0.0025 -snes_rtol 1.0e-10

// for lev in 0 1 2 3 4; do ./advectlayer -snes_fd -al_exactinit -al_noshow -al_dt 0.01 -al_steps 10 -da_refine $lev | grep error; done

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
            f0;   // scale for source term f(x)
  PetscInt  scheme;
} AppCtx;


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
        x       = dx/2 + dx * (PetscReal)j;
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
  return user->v0;
}

// alternative for velocity():
//   const PetscReal Lsqr = user->L * user->L;
//   return (4.0 / Lsqr) * user->v0 * x * (user->L - x);


// without constraint, with this f(x), \int_0^L u(t,x) dt --> - \infty
PetscReal fsource(const PetscReal x, const AppCtx *user) {
  const PetscReal CC = 2.0 * PETSC_PI / user-> L;
  return - (user->f0/5.0) + user->f0 * sin(CC * x);
}


/* for call-back: evaluate residual FF(x) on local process patch */
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info,PetscScalar *u,PetscScalar *FF,
                                 AppCtx *user) {
  PetscErrorCode  ierr;
  PetscInt        j;
  const PetscReal dt = user->dt, dx = user->dx,
                  nu = dt / dx, nu2 = nu / 2.0;
  PetscReal       *uold;

  PetscFunctionBeginUser;
  ierr = DMDAVecGetArray(info->da, user->uold, &uold);CHKERRQ(ierr);
  for (j=info->xs; j<info->xs+info->xm; j++) {
      const PetscReal x = dx/2 + dx * (PetscReal)j;
      switch (user->scheme) {
        case 0 : { // backward-Euler, centered scheme; works even if v=v(x)
          const PetscReal vleft  = velocity(x - dx/2,user),
                          vright = velocity(x + dx/2,user);
          FF[j] = - nu2 * vleft * u[j-1]
                  + (1.0 + nu2 * (vright - vleft)) * u[j]
                  + nu2 * vright * u[j+1]
                  - uold[j] - dt * fsource(x,user);
          break; }
        case 1 : { // this third-order upwind biased implicit thing only works if v(x)=v0>0
          const PetscReal mu = user->v0 * nu / 6.0;
          FF[j] =   mu * u[j-2]
                  - 6.0*mu * u[j-1]
                  + (1.0 + 3.0*mu) * u[j]
                  + 2.0*mu * u[j+1]
                  - uold[j] - dt * fsource(x,user);
          break; }
        case 2 : { // the box scheme; works even if v=v(x)
          const PetscReal vleft  = velocity(x - dx/2,user),
                          vright = velocity(x + dx/2,user);
          FF[j] = u[j-1] + u[j] - uold[j-1] - uold[j]
                  + nu * ( vright * (u[j] + uold[j]) - vleft * (u[j-1] + uold[j-1]) )
                  - 2.0 * dt * fsource(x,user);
          break; }
        default :
          SETERRQ(PETSC_COMM_WORLD,2,"not allowed value of scheme\n");
      }
  }
  ierr = DMDAVecRestoreArray(info->da, user->uold, &uold);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


int main(int argc,char **argv) {
  PetscErrorCode      ierr;
  SNES                snes;
  Vec                 u;
  AppCtx              user;
  PetscInt            NN;
  PetscBool           noshow, useexactinit, o3set, boxset;
  DMDALocalInfo       info;

  PetscInitialize(&argc,&argv,(char*)0,help);

  user.L  = 100.0;
  user.v0 = 100.0;
  user.f0 = 1.0;

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"al_","options to advectlayer","");CHKERRQ(ierr);
  NN = 10;
  ierr = PetscOptionsInt("-steps","number of time steps",
                         NULL,NN,&NN,NULL);CHKERRQ(ierr);
  user.dt = 0.05;
  ierr = PetscOptionsReal("-dt","length of time step",
                          NULL,user.dt,&user.dt,NULL);CHKERRQ(ierr);
  noshow = PETSC_FALSE;
  ierr = PetscOptionsBool("-noshow","do not show solution with X",
                          NULL,noshow,&noshow,NULL);CHKERRQ(ierr);
  useexactinit = PETSC_FALSE;
  ierr = PetscOptionsBool("-exactinit","initialize with exact solution",
                          NULL,useexactinit,&useexactinit,NULL);CHKERRQ(ierr);
  o3set = PETSC_FALSE;
  boxset = PETSC_FALSE;
  ierr = PetscOptionsBool("-thirdorder","use third-order upwind-biased FD scheme",
                          NULL,o3set,&o3set,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-box","use box method FD scheme",
                          NULL,boxset,&boxset,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  if (o3set && boxset) {
    SETERRQ(PETSC_COMM_WORLD,1,"set -thirdorder or -box or neither, but BOTH -thirdorder and -box is conflict\n");
  } else if (o3set)
    user.scheme = 1;
  else if (boxset)
    user.scheme = 2;
  else
    user.scheme = 0;

  ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC,
                      -50,         // override with -da_grid_x or -da_refine
                      1, 2, NULL,  // dof = 1 and stencil width = 2
                      &user.da); CHKERRQ(ierr);
  ierr = DMSetApplicationContext(user.da, &user);CHKERRQ(ierr);

  ierr = DMDAGetLocalInfo(user.da,&info); CHKERRQ(ierr);
  user.dx = user.L / (PetscReal)(info.mx);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "doing %d steps of dt=%f on grid of %d points with spacing %g\n"
                     "[CFL time step is %g]\n",
                     NN,user.dt,info.mx,user.dx,user.dx/user.v0); CHKERRQ(ierr);
  // cell-centered grid
  ierr = DMDASetUniformCoordinates(user.da,user.dx/2,user.L-user.dx/2,
                                   -1.0,-1.0,-1.0,-1.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(user.da,&u);CHKERRQ(ierr);
  ierr = VecSetOptionsPrefix(u,"u_"); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)u,"solution u"); CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(user.da,&user.uold);CHKERRQ(ierr);
  ierr = VecSetOptionsPrefix(user.uold,"uold_"); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)user.uold,"OLD solution uold"); CHKERRQ(ierr);

  if (useexactinit) {
    ierr = SetToExactSolution(user.uold,&user);CHKERRQ(ierr);
  } else {
    ierr = VecSet(user.uold,0.0); CHKERRQ(ierr);
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
    for (n = 0; n < NN; ++n) {
      ierr = PetscPrintf(PETSC_COMM_WORLD, "  time %7g: ", t+user.dt); CHKERRQ(ierr);
      ierr = VecScale(u,0.99); CHKERRQ(ierr); // move u so that line search goes some distance
      ierr = SNESSolve(snes, NULL, u); CHKERRQ(ierr);
      ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
      ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"%3d Newton iterations;   result = %s\n",
                         its,SNESConvergedReasons[reason]);CHKERRQ(ierr);
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
  Vec uexact;
  PetscScalar errnorm;
  ierr = DMCreateGlobalVector(user.da,&uexact);CHKERRQ(ierr);
  ierr = SetToExactSolution(uexact,&user);CHKERRQ(ierr);
  ierr = VecAXPY(u,-1.0,uexact); CHKERRQ(ierr);    // u <- u + (-1.0) uxact
  ierr = VecNorm(u,NORM_INFINITY,&errnorm); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
             "on dx=%.4e grid with N=%4d steps of dt=%.4e:  error |u-uexact|_inf = %g\n",
             user.dx,NN,user.dt,errnorm); CHKERRQ(ierr);

  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&user.uold);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = DMDestroy(&user.da);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}

