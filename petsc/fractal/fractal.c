static char help[] =
"Structured-grid *decoupled* Newton solver test problem.\n\n";

// FIXME

/*
the reason I wrote this odd-ball (and stub) code is that I am thinking about the sia-fve project and its convergence difficulties.  question is whether decoupled residuals like
   f_{ij}(u) = (u_{ij} + b(x_i,y_j))^5 + u_{ij}
all with same initial iterate  u_{ij}^{(0)} = 0,  but having different functions because of highly-variable  b(x,y),  have convergence difficulties anything like the complex newton fractals that come from solving  f(z) = 0  for grid of z_0 = (x_i,y_j)
*/

//  ./fractal -snes_monitor -snes_rtol 1.0e-16 -da_refine 4 -snes_converged_reason -snes_stol 1.0e-16 -snes_atol 0.0 -snes_fd_color


#include <petsc.h>

typedef struct {
  PetscViewer viewer;
} Ctx;

PetscErrorCode FormFunctionLocal(DMDALocalInfo *info, double **au,
                                 double **FF, void *user) {
    const double hx = 10.0/(info->mx-1),  hy = 10.0/(info->my-1);
    int          i, j;
    double       x, y, b;

    for (j = info->ys; j < info->ys + info->ym; j++) {
        y = j * hy;
        for (i = info->xs; i < info->xs + info->xm; i++) {
            x = i * hx;
            b = x + 10.0 * sin(y) * sin(5.0 * x);
            FF[j][i] = pow(au[j][i] + b,5.0) + au[j][i];
        }
    }
    return 0;
}

/*
PetscErrorCode FormJacobianLocal(DMDALocalInfo *info, PetscScalar **au,
                                 Mat J, Mat Jpre, void *user) {
    PetscErrorCode  ierr;
    int          i,j,;
    double       v;
    MatStencil   col,row;

    for (j = info->ys; j < info->ys + info->ym; j++) {
        y = j * hy;
        for (i = info->xs; i < info->xs + info->xm; i++) {
            x = i * hx;
            b = x + 10.0 * sin(y) * sin(5.0 * x);
            row =
            col =
            v = 
            ierr = MatSetValuesStencil(Jpre,1,&row,1,&col,&v,INSERT_VALUES); CHKERRQ(ierr);
        }
    }

    ierr = MatAssemblyBegin(Jpre,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jpre,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    if (J != Jpre) {
        ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    }
    return 0;
}
*/

PetscErrorCode CountMonitor(SNES snes, PetscInt its, double fnorm, void *ctx)
{
  PetscErrorCode ierr;
  Ctx            *user = (Ctx*) ctx;
  Vec            x;

  //ierr = PetscPrintf(PETSC_COMM_WORLD,
  //                   "iter %D:  residual norm = %g\n",its,fnorm); CHKERRQ(ierr);
  ierr = SNESGetSolutionUpdate(snes,&x);CHKERRQ(ierr);
  ierr = VecView(x,user->viewer);CHKERRQ(ierr);
  return 0;
}

int main(int argc,char **argv) {
  PetscErrorCode ierr;
  SNES           snes;
  Vec            u;
  DM             da;
  DMDALocalInfo  info;
  Ctx            user;

  PetscInitialize(&argc,&argv,NULL,help);

  ierr = DMDACreate2d(PETSC_COMM_WORLD,
               DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,
               -9,-9,PETSC_DECIDE,PETSC_DECIDE,1,0,NULL,NULL,&da); CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(da,&u);CHKERRQ(ierr);

  ierr = VecSet(u,0.0); CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetDM(snes,da); CHKERRQ(ierr);
  ierr = DMDASNESSetFunctionLocal(da,INSERT_VALUES,
             (DMDASNESFunction)FormFunctionLocal,NULL); CHKERRQ(ierr);
//  ierr = DMDASNESSetJacobianLocal(da,
//             (DMDASNESJacobian)FormJacobianLocal,&user); CHKERRQ(ierr);
  ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,0,0,0,0,400,400,&user.viewer);CHKERRQ(ierr);
//  ierr = SNESMonitorSet(snes,CountMonitor,&user,NULL); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  ierr = SNESSolve(snes,NULL,u); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"done on %d x %d grid\n",info.mx,info.my); CHKERRQ(ierr);

  VecDestroy(&u);
  SNESDestroy(&snes);  DMDestroy(&da);
  PetscViewerDestroy(&user.viewer);
  PetscFinalize();
  return 0;
}

