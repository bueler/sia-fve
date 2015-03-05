/* (C) 2015 Ed Bueler */

#include <petscdmda.h>
#include <petscvec.h>

#include "mahaffyctx.h"
#include "exactsia.h"


PetscErrorCode GetCoordArray(DMDACoor2d ***acoord, const AppCtx *user) {
  PetscErrorCode ierr;
  DM            coordDA;
  Vec           coordinates;

  PetscFunctionBeginUser;
  ierr = DMGetCoordinateDM(user->da, &coordDA); CHKERRQ(ierr);
  ierr = DMGetCoordinates(user->da, &coordinates); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(coordDA, coordinates, acoord); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode RestoreCoordArray(DMDACoor2d ***acoord, const AppCtx *user) {
  PetscErrorCode ierr;
  DM            coordDA;
  Vec           coordinates;

  PetscFunctionBeginUser;
  ierr = DMGetCoordinateDM(user->da, &coordDA); CHKERRQ(ierr);
  ierr = DMGetCoordinates(user->da, &coordinates); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(coordDA, coordinates, acoord); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#define domeL  750.0e3 // radius of exact ice sheet
#define domeH0 3600.0  // center thickness of exact ice sheet

PetscReal radialcoord(const DMDACoor2d c) {
  PetscReal r;
  r = PetscSqrtReal(c.x * c.x + c.y * c.y);
  if (r < 0.01)
      r = 0.01;  // avoid r=0 singularity in dome formulas
  return r;
}

PetscErrorCode SetToDomeCMB(Vec m, const AppCtx *user) {
  PetscErrorCode ierr;
  const PetscReal L  = domeL,
                  H0 = domeH0,
                  n  = user->n,
                  pp = 1.0 / n,
                  CC = user->Gamma * PetscPowReal(H0,2.0*n+2.0)
                          / PetscPowReal(2.0 * L * (1.0-1.0/n),n);
  DMDALocalInfo info;
  DMDACoor2d    **coords;
  PetscReal     **am, r, s, tmp1, tmp2;
  PetscInt      j, k;

  PetscFunctionBeginUser;
  ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
  ierr = GetCoordArray(&coords, user); CHKERRQ(ierr);
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
  ierr = RestoreCoordArray(&coords, user); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode SetToDomeExactThickness(Vec H, const AppCtx *user) {
  PetscErrorCode ierr;
  const PetscReal L  = domeL,
                  H0 = domeH0,
                  n  = user->n,
                  mm = 1.0 + 1.0 / n,
                  qq = n / (2.0 * n + 2.0),
                  CC = H0 / PetscPowReal(1.0 - 1.0 / n,qq);
  DMDALocalInfo info;
  DMDACoor2d    **coords;
  PetscReal     **aH, r, s, tmp;
  PetscInt      j, k;

  PetscFunctionBeginUser;
  ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
  ierr = GetCoordArray(&coords, user); CHKERRQ(ierr);
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
  ierr = RestoreCoordArray(&coords, user); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// following constants are from page 236 of Jarosch et al (2013), below eqn (62)
#define secpera   3.1556926e7
#define bedstepxm 20000.0
#define bedstepxs 7000.0
#define bedstepb0 500.0
#define bedstepm0 2.0/secpera

// formula (45) in Jarosch et al (2013)
PetscErrorCode SetToBedStepBed(Vec b, const AppCtx *user) {
  PetscErrorCode ierr;
  DMDALocalInfo info;
  DMDACoor2d    **coords;
  PetscReal     **ab, absx;
  PetscInt      j, k;

  PetscFunctionBeginUser;
  ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
  ierr = GetCoordArray(&coords, user); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, b, &ab);CHKERRQ(ierr);
  for (k=info.ys; k<info.ys+info.ym; k++) {
      for (j=info.xs; j<info.xs+info.xm; j++) {
          absx = PetscAbsReal(coords[k][j].x);
          if (absx < bedstepxs)
              ab[k][j] = bedstepb0;
          else
              ab[k][j] = 0.0;
      }
  }
  ierr = DMDAVecRestoreArray(user->da, b, &ab);CHKERRQ(ierr);
  ierr = RestoreCoordArray(&coords, user); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


// formula (54) in Jarosch et al (2013)
PetscErrorCode SetToBedStepCMB(Vec m, const AppCtx *user) {
  PetscErrorCode ierr;
  const PetscReal n  = user->n,
                  xm = bedstepxm,
                  C  = n * bedstepm0 / PetscPowReal(xm,2.0*n-1.0);
  DMDALocalInfo info;
  DMDACoor2d    **coords;
  PetscReal     **am, absx;
  PetscInt      j, k;

  PetscFunctionBeginUser;
  ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
  ierr = GetCoordArray(&coords, user); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, m, &am);CHKERRQ(ierr);
  for (k=info.ys; k<info.ys+info.ym; k++) {
      for (j=info.xs; j<info.xs+info.xm; j++) {
          absx = PetscAbsReal(coords[k][j].x);
          am[k][j] = C * PetscPowReal(absx * PetscAbsReal(xm - absx),n-1) * (xm - 2.0 * absx);
          am[k][j] = PetscMax(am[k][j], - bedstepm0);  // limit below at -2 m/a
      }
  }
  ierr = DMDAVecRestoreArray(user->da, m, &am);CHKERRQ(ierr);
  ierr = RestoreCoordArray(&coords, user); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


// formulas (56)--(59) in Jarosch et al (2013)
// Note  "h_{s+}(x)" in (58) is really  h_{s+} = lim_{x \to x_s^+} h(x), so it does
//   not depend on x.  Similar comments apply to (59).
PetscErrorCode SetToBedStepExactThickness(Vec Hexact, const AppCtx *user) {
  PetscErrorCode ierr;
  const PetscReal xm = bedstepxm,
                  xs = bedstepxs,
                  b0 = bedstepb0,
                  m0 = bedstepm0,
                  n  = user->n,
                  A  = user->A,
                  rg = user->rho * user->g,
                  ninv    = 1.0 / n,
                  p       = n / (2.0*n+2.0),
                  pinv    = 1.0 / p,
                  q       = (2.0*n-1.0) / n,
                  denom   = 6.0 * n * rg * PetscPowReal(2.0*A,ninv) * PetscPowReal(xm,q),
                  Ch      = (2.0*n+2.0) * PetscPowReal((n+2.0)*m0,ninv) / denom,
                  hsplus  = PetscPowReal(Ch * (xm + 2.0*xs) * (xm - xs) * (xm - xs), p),
                  hsminus = PetscMax(hsplus - b0, 0.0),
                  Ca      = PetscPowReal(hsminus, pinv) - PetscPowReal(hsplus, pinv);
  DMDALocalInfo info;
  DMDACoor2d    **coords;
  PetscReal     **aHex, absx, tmp;
  PetscInt      j, k;

  PetscFunctionBeginUser;
  ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
  ierr = GetCoordArray(&coords, user); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, Hexact, &aHex);CHKERRQ(ierr);
  for (k=info.ys; k<info.ys+info.ym; k++) {
      for (j=info.xs; j<info.xs+info.xm; j++) {
          absx = PetscAbsReal(coords[k][j].x);
          if (absx >= xm)
              aHex[k][j] = 0.0;
          else {
              tmp = Ch * (xm + 2.0*absx) * (xm - absx) * (xm - absx);
              if (absx < xs)
                  aHex[k][j] = PetscPowReal(Ca + tmp,p);
              else
                  aHex[k][j] = PetscPowReal(tmp,p);
          }
      }
  }
  ierr = DMDAVecRestoreArray(user->da, Hexact, &aHex);CHKERRQ(ierr);
  ierr = RestoreCoordArray(&coords, user); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

