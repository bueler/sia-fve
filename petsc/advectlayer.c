static const char help[] = "Solves advecting-layer problem in 1d:\n"
"    u_t + div q = f\n"
"where\n"
"    q = v(x) u\n"
"on 0 < x < L, with periodic boundary conditions, subject to constraint\n"
"    u >= 0.\n"
"Implicit (backward Euler) centered finite difference method and\n"
"SNESVI.\n\n";

//FIXME: add Jacobian

//FIXME: use 3rd-order upwinding with flux-limiter

//FIXME: for 3rd-order, Jacobian ignors flux limiter?

//FIXME: relate to event detection in \infty dimensions

//   ./advectlayer -help |grep al_

//   ./advectlayer -snes_fd -draw_pause 0.5
//   ./advectlayer -snes_mf -draw_pause 0.5

//   ./advectlayer -snes_fd -al_noshow
//   ./advectlayer -snes_fd -al_steps 100
//   ./advectlayer -snes_fd -al_dt 0.1

//   ./advectlayer -snes_fd -snes_type vinewtonssls
//   ./advectlayer -snes_fd -snes_vi_monitor

//   ./advectlayer -snes_fd -al_dt 0.01 -al_steps 1000 -da_refine 2

//   ./advectlayer -snes_fd -al_steps 1000 -da_refine 3 -al_dt 0.0025 -snes_rtol 1.0e-10

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
} AppCtx;


//  for call-back: tell SNESVI (variational inequality) that we want  0.0 <= u < +infinity
PetscErrorCode FormBounds(SNES snes, Vec Xl, Vec Xu)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = VecSet(Xl,0.0); CHKERRQ(ierr);
  ierr = VecSet(Xu,PETSC_INFINITY); CHKERRQ(ierr);
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
  const PetscReal CC = 4.0 * PETSC_PI;
  return -(user->f0/5.0) + user->f0 * sin(CC * x / user->L);
}


/* for call-back: evaluate nonlinear function FF(x) on local process patch */
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info,PetscScalar *u,PetscScalar *FF,
                                 AppCtx *user) {
  PetscErrorCode  ierr;
  PetscInt        j;
  const PetscReal dt = user->dt, dx = user->dx,
                  nu = dt / dx, nu2 = nu / 2.0;
  PetscReal       *uold, x, vleft, vright;

  PetscFunctionBeginUser;
  ierr = DMDAVecGetArray(info->da, user->uold, &uold);CHKERRQ(ierr);
  for (j=info->xs; j<info->xs+info->xm; j++) {
      x      = dx/2 + dx * (PetscReal)j;
      vleft  = velocity(x - dx/2,user);
      vright = velocity(x + dx/2,user);
      // compute residual
      FF[j] = - nu2 * vleft * u[j-1]
              + (1.0 + nu2 * (vright - vleft)) * u[j]
              + nu2 * vright * u[j+1]
              - uold[j] - dt * fsource(x,user);
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
  PetscBool           noshow;
  DM                  da;
  DMDALocalInfo       info;

  PetscInitialize(&argc,&argv,(char*)0,help);

  user.L  = 10000.0;
  user.v0 = 10000.0;
  user.f0 = 100.0;

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
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC,
                      -50,         // override with -da_grid_x or -da_refine
                      1, 1, NULL,  // dof = 1 and stencil width = 1
                      &da); CHKERRQ(ierr);
  ierr = DMSetApplicationContext(da, &user);CHKERRQ(ierr);

  ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
  user.dx = user.L / (PetscReal)(info.mx);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "doing %d steps of dt=%f on grid of %d points with spacing %g\n"
                     "[CFL time step is %g]\n",
                     NN,user.dt,info.mx,user.dx,user.dx/user.v0); CHKERRQ(ierr);
  // cell-centered grid
  ierr = DMDASetUniformCoordinates(da,user.dx/2,user.L-user.dx/2,
                                   -1.0,-1.0,-1.0,-1.0);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(da,&u);CHKERRQ(ierr);
  ierr = VecSetOptionsPrefix(u,"u_"); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)u,"solution u"); CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(da,&user.uold);CHKERRQ(ierr);
  ierr = VecSetOptionsPrefix(user.uold,"uold_"); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)user.uold,"OLD solution uold"); CHKERRQ(ierr);

  ierr = VecSet(user.uold,0.0);CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetDM(snes,da);CHKERRQ(ierr);
  ierr = SNESSetType(snes,SNESVINEWTONRSLS);CHKERRQ(ierr);
  ierr = SNESVISetComputeVariableBounds(snes,&FormBounds);CHKERRQ(ierr);

  ierr = DMDASNESSetFunctionLocal(da,INSERT_VALUES,
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
      ierr = SNESSolve(snes, NULL, u); CHKERRQ(ierr);
      ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
      ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"%3d Newton iterations;   result = %s\n",
                         its,SNESConvergedReasons[reason]);CHKERRQ(ierr);
      if (!noshow) {
        ierr = VecView(u, PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);
      }
      ierr = VecCopy(u, user.uold); CHKERRQ(ierr);
      t += user.dt;
    }
  }

  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&user.uold);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}

