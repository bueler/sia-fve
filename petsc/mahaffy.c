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

//   ./mahaffy -help |grep mah_

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
  PetscBool star,   // use Mahaffy* instead of Mahaffy
            exactinit;// initialize with exact solution instead of zero
} AppCtx;


extern PetscErrorCode FormBounds(SNES,Vec,Vec);
extern PetscErrorCode SetToExactThickness(Vec,const AppCtx*);
//extern PetscErrorCode SetToSMB(Vec,const AppCtx*);
extern PetscErrorCode ProcessOptions(AppCtx*);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*,PetscScalar**,PetscScalar**,AppCtx*);
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
  
  user.star   = PETSC_FALSE;
  user.exactinit = PETSC_FALSE;

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

  if (user.exactinit == PETSC_TRUE) {
      ierr = SetToExactThickness(H,&user);CHKERRQ(ierr);
  } else {
      ierr = VecSet(H,0.0); CHKERRQ(ierr);
  }

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


PetscReal radialcoord(const DMDACoor2d c) {
  PetscReal r;
  r = PetscSqrtReal(c.x * c.x + c.y * c.y);
  if (r < 0.01)
      r = 0.01;  // avoid r=0 singularity
  return r;
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


PetscErrorCode SetToExactThickness(Vec H, const AppCtx *user) {
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
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


/*
PetscReal SetToSMB(Vec m, const AppCtx *user) {
  FIXME
  return FIXME;
}
*/

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
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info,PetscScalar **H,PetscScalar **FF,
                                 AppCtx *user) {
  //PetscErrorCode  ierr;
  //const PetscReal dx = user->dx;
  PetscInt        j, k;

  PetscFunctionBeginUser;
  for (k=info->ys; k<info->ys+info->ym; k++) {
      for (j=info->xs; j<info->xs+info->xm; j++) {
          FF[k][j] = 0.0 * H[k][j];  //FIXME
      }
  }
  PetscFunctionReturn(0);
}


PetscErrorCode ProcessOptions(AppCtx *user) {
  PetscErrorCode ierr;

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"mah_","options to mahaffy","");CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-exactinit", "initialize with exact solution",
      NULL,user->exactinit,&user->exactinit,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}



