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

//   ./mahaffy -mah_exactinit -help |grep mah_

// views:
//   ./mahaffy -mah_draw -draw_pause 1 -mah_exactinit -da_grid_x 45 -da_grid_y 45
//   ./mahaffy -mah_draw -draw_pause 1 -mah_exactinit

// diverge:
//   ./mahaffy -mah_draw -draw_pause 1
//   ./mahaffy -mah_draw -draw_pause 1 -mah_exactinit -da_grid_x 36 -da_grid_y 36

#include <math.h>
#include <petscdmda.h>
#include <petscsnes.h>

typedef struct {
  DM        da;
  Vec       b,      // the bed elevation
            m;      // the (steady) surface mass balance
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
            exactinit,// initialize with exact solution instead of zero
            draw;   // use X to draw solution H, smb m, and error H-Hexact
} AppCtx;


extern PetscErrorCode FormBounds(SNES,Vec,Vec);
extern PetscErrorCode SetToExactThickness(Vec,const AppCtx*);
extern PetscErrorCode SetToSMB(Vec,const AppCtx*);
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
  user.A      = 1.0e-16/user.secpera; // = 3.17e-24  1/(Pa^3 s); EISMINT I value
  user.Gamma  = 2.0 * PetscPowReal(user.rho*user.g,user.n) * user.A / (user.n+2.0);

  user.exactL = 750.0e3;    // m
  user.exactH0= 3600.0;     // m
  
  user.star      = PETSC_FALSE;
  user.exactinit = PETSC_FALSE;
  user.draw      = PETSC_FALSE;

  ierr = ProcessOptions(&user); CHKERRQ(ierr);

  ierr = DMDACreate2d(PETSC_COMM_WORLD,
                      DM_BOUNDARY_PERIODIC,DM_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX,
                      -9,-9,PETSC_DECIDE,PETSC_DECIDE,  // default 9x9 grid has dx=200km
                      1, 1,                             // dof=1,  stencilwidth=1
                      NULL,NULL,&user.da);
  ierr = DMSetApplicationContext(user.da, &user);CHKERRQ(ierr);

  ierr = DMDAGetLocalInfo(user.da,&info); CHKERRQ(ierr);
  if (info.mx != info.my)  {
      SETERRQ(PETSC_COMM_WORLD,1,"mx != my ERROR: grid must have equal spacing in x,y ...\n");
  }

  user.dx = 2 * user.L / (PetscReal)(info.mx);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
                     "solving on %d x %d grid spacing dx=%g ...\n",
                     info.mx,info.my,user.dx); CHKERRQ(ierr);

  // cell-centered grid
  ierr = DMDASetUniformCoordinates(user.da,
             -user.L+user.dx/2,user.L-user.dx/2,
             -user.L+user.dx/2,user.L-user.dx/2,
             0.0,1.0); CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(user.da,&H);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)H,"thickness solution H"); CHKERRQ(ierr);
  if (user.exactinit == PETSC_TRUE) {
      ierr = SetToExactThickness(H,&user);CHKERRQ(ierr);
  } else {
      ierr = VecSet(H,0.0); CHKERRQ(ierr);
  }

  ierr = VecDuplicate(H,&user.b); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)(user.b),"bed elevation b"); CHKERRQ(ierr);
  ierr = VecSet(user.b,0.0); CHKERRQ(ierr);

  ierr = VecDuplicate(H,&user.m); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)(user.m),"surface mass balance m"); CHKERRQ(ierr);
  ierr = SetToSMB(user.m,&user);CHKERRQ(ierr);

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

  if (user.draw == PETSC_TRUE) {
      ierr = VecView(H, PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);
      ierr = VecView(user.m, PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);
  }

  PetscScalar enorminf,enorm1;
  Vec         Hexact;
  ierr = VecDuplicate(H,&Hexact); CHKERRQ(ierr);
  ierr = SetToExactThickness(Hexact,&user); CHKERRQ(ierr);
  ierr = VecAXPY(H,-1.0,Hexact); CHKERRQ(ierr);    // H < - H + (-1.0) Hexact
  ierr = VecNorm(H,NORM_INFINITY,&enorminf); CHKERRQ(ierr);
  ierr = VecNorm(H,NORM_1,&enorm1); CHKERRQ(ierr);
  enorm1 /= info.mx * info.my;
  ierr = PetscPrintf(PETSC_COMM_WORLD,
             "errors:  max |u-uexact| = %g,  av |u-uexact| = %g\n",enorminf,enorm1); CHKERRQ(ierr);
  ierr = VecDestroy(&Hexact);CHKERRQ(ierr);

  if (user.draw == PETSC_TRUE) {
      ierr = VecView(H, PETSC_VIEWER_DRAW_WORLD); CHKERRQ(ierr);
  }

  ierr = VecDestroy(&user.m);CHKERRQ(ierr);
  ierr = VecDestroy(&user.b);CHKERRQ(ierr);
  ierr = VecDestroy(&H);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
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


PetscErrorCode SetToSMB(Vec m, const AppCtx *user) {
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
      }
  }
  ierr = DMDAVecRestoreArray(user->da, m, &am);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


// single formula for SIA
PetscReal getD(const PetscReal H, const PetscReal sx, const PetscReal sy, const AppCtx *user) {
    PetscReal n = user->n;
    return user->Gamma * PetscPowReal(H,n+2.0) * PetscPowReal(sx*sx + sy*sy,(n-1.0)/2);
}

typedef struct {
    PetscReal x,y;
} Flux;

typedef struct {
    PetscReal x,y;
} Grad;


/* evaluate the flux using Q1 approximations to H(x,y) and b(x,y) on element */
PetscErrorCode fluxatpt(PetscInt j, PetscInt k,         // (j,k) is the element (by lower-left corner)
                        PetscReal locx, PetscReal locy, // = (x,y) coords in element
                        PetscReal **H, PetscReal **b,   // H[k][j] and b[k][j] are node values
                        const AppCtx *user, Flux *q, Grad *grads) {
  const PetscReal dx = user->dx, dy = dx,
                  x[4]  = {1.0 - locx / dx, locx / dx,       1.0 - locx / dx, locx / dx},
                  gx[4] = {- 1.0 / dx,      1.0 / dx,        - 1.0 / dx,      1.0 / dx},
                  y[4]  = {1.0 - locy / dy, 1.0 - locy / dy, locy / dy,       locy / dy},
                  gy[4] = {- 1.0 / dy,      - 1.0 / dy,      1.0 / dy,        1.0 / dy};
  PetscReal HH, sx, sy, DD;
  if ((q == NULL) || (H == NULL) || (b == NULL)) {
      SETERRQ(PETSC_COMM_WORLD,1,"ERROR: illegal NULL ptr in fluxatpt() ...\n");
  }
  HH = x[0]*y[0]*H[k][j] + x[1]*y[1]*H[k][j+1] + x[2]*y[2]*H[k+1][j] + x[3]*y[3]*H[k+1][j+1];
  sx = gx[0]*y[0]*(H[k][j] + b[k][j])
                         + gx[1]*y[1]*(H[k][j+1] + b[k][j+1])
                                               + gx[2]*y[2]*(H[k+1][j] + b[k+1][j])
                                                                     + gx[3]*y[3]*(H[k+1][j+1] + b[k+1][j+1]);
  sy = x[0]*gy[0]*(H[k][j] + b[k][j])
                         + x[1]*gy[1]*(H[k][j+1] + b[k][j+1])
                                               + x[2]*gy[2]*(H[k+1][j] + b[k+1][j])
                                                                     + x[3]*gy[3]*(H[k+1][j+1] + b[k+1][j+1]);
  DD = getD(HH,sx,sy,user);
  q->x = - DD * sx;
  q->y = - DD * sy;
  if (grads != NULL) {
      grads->x = sx;
      grads->y = sy;
  }
  PetscFunctionReturn(0);
}


/* for call-back: evaluate residual FF(x) on local process patch */
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info,PetscScalar **H,PetscScalar **FF,
                                 AppCtx *user) {
  PetscErrorCode  ierr;
  const PetscReal dx = user->dx, dy = dx;
  PetscInt        j, k;
  PetscReal       **am, **ab;
  Flux            qeu, qed, qnl, qnr, qwu, qwd, qsl, qsr;
  Vec             bloc;

  PetscFunctionBeginUser;
  ierr = DMCreateLocalVector(user->da,&bloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(user->da,user->b,INSERT_VALUES,bloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(user->da,user->b,INSERT_VALUES,bloc);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(user->da, bloc, &ab);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->m, &am);CHKERRQ(ierr);
  for (k=info->ys; k<info->ys+info->ym; k++) {
      for (j=info->xs; j<info->xs+info->xm; j++) {
          ierr = fluxatpt(j,   k,       dx/2.0,     dy/4.0, H,ab,user,&qeu,NULL); CHKERRQ(ierr);
          ierr = fluxatpt(j,   k-1,     dx/2.0, 3.0*dy/4.0, H,ab,user,&qed,NULL); CHKERRQ(ierr);
          ierr = fluxatpt(j-1, k,   3.0*dx/4.0,     dy/2.0, H,ab,user,&qnl,NULL); CHKERRQ(ierr);
          ierr = fluxatpt(j,   k,       dx/4.0,     dy/2.0, H,ab,user,&qnr,NULL); CHKERRQ(ierr);
          ierr = fluxatpt(j-1, k,       dx/2.0,     dy/4.0, H,ab,user,&qwu,NULL); CHKERRQ(ierr);
          ierr = fluxatpt(j-1, k-1,     dx/2.0, 3.0*dy/4.0, H,ab,user,&qwd,NULL); CHKERRQ(ierr);
          ierr = fluxatpt(j-1, k-1, 3.0*dx/4.0,     dy/2.0, H,ab,user,&qsl,NULL); CHKERRQ(ierr);
          ierr = fluxatpt(j,   k-1,     dx/4.0,     dy/2.0, H,ab,user,&qsr,NULL); CHKERRQ(ierr);
          FF[k][j] =  (qeu.x + qed.x) * dy/2.0 + (qnl.y + qnr.y) * dx/2.0
                    - (qwu.x + qwd.x) * dy/2.0 - (qsl.y + qsr.y) * dx/2.0
                    - am[k][j] * dx * dy;
      }
  }
  ierr = DMDAVecRestoreArray(user->da, user->m, &am);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, bloc, &ab);CHKERRQ(ierr);

  ierr = VecDestroy(&bloc); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode ProcessOptions(AppCtx *user) {
  PetscErrorCode ierr;

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"mah_","options to mahaffy","");CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-exactinit", "initialize with exact solution",
      NULL,user->exactinit,&user->exactinit,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool(
      "-draw", "use X to draw H, m, H-Hexact",
      NULL,user->draw,&user->draw,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}



