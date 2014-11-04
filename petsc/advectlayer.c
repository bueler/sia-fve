static const char help[] = "Solves advecting-layer problem in 1d:\n"
"  u_t + div q = f  where  q = v(x) u\n"
"on 0 < x < L but subject to constraint\n"
"  u >= 0.\n";

//  ./advectlayer -help |grep ad_

//  ./advectlayer -ad_steps 100
//  ./advectlayer -ad_dt 0.1

//  ./advectlayer -snes_type vinewtonssls
//  ./advectlayer -snes_vi_monitor

#include <math.h>
#include <petscdmda.h>
#include <petscsnes.h>


typedef struct {
  DM        da;
  Vec       uold;
  PetscReal dt,
            dx,
            L,
            v0,
            f0,
            q0,
            qL;
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


/* for call-back: evaluate nonlinear function FF(x) on local process patch */
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info,PetscScalar *u,PetscScalar *FF,AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       j;
  PetscReal      *uold, 
                 x, xleft, xright, vleft, vright, f,
                 nu   = user->dt / user->dx,
                 nu2  = nu / 2.0,
                 pi   = PETSC_PI,
                 Lsqr = user->L * user->L,
                 L5   = user->L / 5.0;

  PetscFunctionBeginUser;
  ierr = DMDAVecGetArray(info->da, user->uold, &uold);CHKERRQ(ierr);
  for (j=info->xs; j<info->xs+info->xm; j++) {
      x      = user->dx * ((PetscReal)j + 1.0/2.0);
      xleft  = x - user->dx/2;
      xright = x + user->dx/2;
      vleft  = (4.0 / Lsqr) * user->v0 * (user->L - xleft);
      vright = (4.0 / Lsqr) * user->v0 * (user->L - xright);
      //vleft  = user->v0;
      //vright = user->v0;
      if (x <= L5 || x >= 4.0*L5)
          f = 0.0;
      else
          f = user->f0 * sin(3.0 * pi * (x - L5) / (3.0*L5));
      if (j == 0) {
          FF[j] = - nu * user->q0
                  + (1.0 + nu2 * vright) * u[j]
                  + nu2 * vright * u[j+1];
      } else if (j == info->mx-1) {
          FF[j] = - nu2 * vleft * u[j-1]
                  + (1.0 - nu2 * vleft) * u[j]
                  + nu * user->qL;
      } else {
          FF[j] = - nu2 * vleft * u[j-1]
                  + (1.0 + nu2 * (vright - vleft)) * u[j]
                  + nu2 * vright * u[j+1];
      }
      FF[j] -= uold[j] + user->dt * f;
  }
  ierr = DMDAVecRestoreArray(info->da, user->uold, &uold);CHKERRQ(ierr);

  //FIXME ierr = PetscLogFlops(10.0*info->ym*info->xm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


int main(int argc,char **argv)
{
  PetscErrorCode      ierr;
  SNES                snes;
  Vec                 u;
  AppCtx              user;
  PetscReal           t;
  PetscInt            its, n, NN;
  SNESConvergedReason reason;
  DM                  da;
  DMDALocalInfo       info;

  PetscInitialize(&argc,&argv,(char*)0,help);

  user.L  = 10000.0;
  user.v0 = 10000.0;
  user.f0 = 100.0;
  user.q0 = 0.0;
  user.qL = 0.0;

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"ad_","options to advectlayer","");CHKERRQ(ierr);
  NN = 10;
  ierr = PetscOptionsInt("-steps","number of time steps",NULL,NN,&NN,NULL);CHKERRQ(ierr);
  user.dt = 30.0;
  ierr = PetscOptionsReal("-dt","length of time step",NULL,user.dt,&user.dt,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE,
                      -100,        // override with -da_grid_x or -da_refine
                      1, 1, NULL,  // dof = 1 and stencil width = 1
                      &da); CHKERRQ(ierr);
  ierr = DMSetApplicationContext(da, &user);CHKERRQ(ierr);

  ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
  user.dx = user.L / (PetscReal)(info.mx);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "using grid of %d points with spacing %g ...\n",
                     info.mx,user.dx); CHKERRQ(ierr);
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

  /* Time Loop */
  t = 0.0;
  ierr = VecCopy(user.uold,u);CHKERRQ(ierr);
  for (n = 0; n < NN; ++n) {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "  time %7g: ", t); CHKERRQ(ierr);
    ierr = VecView(u, PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);

    ierr = SNESSolve(snes, NULL, u); CHKERRQ(ierr);
    ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
    ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"%3d Newton iterations;   result = %s\n",
                       its,SNESConvergedReasons[reason]);CHKERRQ(ierr);

    ierr = VecCopy(u, user.uold); CHKERRQ(ierr);
    ierr = VecView(u, PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);
    t += user.dt;
  }

  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&user.uold);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}

