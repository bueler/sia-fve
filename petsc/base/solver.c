/* (C) 2015 Ed Bueler */

#include <petscsnes.h>

#include "q1op.h"
#include "sia.h"
#include "io.h"  // for myPrintf()
#include "solver.h"

/* Loop over locally-owned elements, including ghosts, checking
   nonnegativity of thickness.  Stops with error if not.  */
PetscErrorCode checkadmissible(DMDALocalInfo *info, PetscScalar **H) {
  PetscInt        j, k;
  PetscFunctionBeginUser;
  for (k = info->ys-1; k < info->ys + info->ym + 1; k++) {
      for (j = info->xs-1; j < info->xs + info->xm + 1; j++) {
          if (H[k][j] < 0.0)
              SETERRQ3(PETSC_COMM_WORLD,1,"ERROR: inadmissible H[%d][%d] = %.3e < 0\n",k,j,H[k][j]);
      }
  }
  PetscFunctionReturn(0);
}


// averages two gradients; only used in true Mahaffy
Grad gradav(Grad g1, Grad g2) {
  Grad gav;
  gav.x = (g1.x + g2.x) / 2.0;
  gav.y = (g1.y + g2.y) / 2.0;
  return gav;
}


/* For call-back by SNES using DMDA info.

Evaluates residual FF on local process patch:
   FF_{j,k} = \int_{\partial V_{j,k}} \mathbf{q} \cdot \mathbf{n} - m_{j,k} \Delta x \Delta y
where V_{j,k} is the control volume centered at (x_j,y_k).

Regarding indexing the location along the boundary of the control volume where
flux is evaluated, this shows four elements and one control volume centered
at (x_j,y_k).  The boundary of the control volume has 8 points, numbered s=0,...,7:
   -------------------
  |         |         |
  |    ..2..|..1..    |
  |   3:    |    :0   |
k |--------- ---------|
  |   4:    |    :7   |
  |    ..5..|..6..    |
  |         |         |
   -------------------
            j

Regarding flux-component indexing on the element indexed by (j,k) node, as shown,
the value  aq[k][j][c]  for c=0,1,2,3, is an x-component at "*" and a y-component
at "%":
   -------------------
  |         :         |
  |         *2        |
  |    3    :    1    |
  |....%.... ....%....|
  |         :         |
  |         *0        |
  |         :         |
  @-------------------
(j,k)
*/
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info, PetscScalar **aH, PetscScalar **FF, AppCtx *user) {
  PetscErrorCode  ierr;
  const PetscReal dx = user->dx, dy = dx;
  FLUXINTCOEFFS
  const PetscBool upwind = (user->lambda > 0.0);
  const PetscReal upmin = (1.0 - user->lambda) * 0.5,
                  upmax = (1.0 + user->lambda) * 0.5;
  PetscInt        j, k;
  PetscReal       **am, **ab, **aHprev, **aDnodemax, **aWmagnodemax,
                  ***aqquad,  ***aDquad, ***aWquad;
  Vec             qloc, Dloc, Wloc;
  DiagnosticScheme *ds = &(user->ds);

  PetscFunctionBeginUser;
  if (user->checkadmissible) {
      ierr = checkadmissible(info,aH); CHKERRQ(ierr);
  }

  if (!user->nodiag) {
      ds->avD      = 0.0;
      ds->avDcount = 0;
      ds->maxD     = 0.0;
      ierr = DMGetLocalVector(user->quadda,&Dloc);CHKERRQ(ierr);
      ierr = DMGetLocalVector(user->quadda,&Wloc);CHKERRQ(ierr);
      ierr = DMDAVecGetArray(user->da, ds->Dnodemax, &aDnodemax);CHKERRQ(ierr);
      ierr = DMDAVecGetArray(user->da, ds->Wmagnodemax, &aWmagnodemax);CHKERRQ(ierr);
      ierr = DMDAVecGetArrayDOF(user->quadda, Dloc, &aDquad);CHKERRQ(ierr);
      ierr = DMDAVecGetArrayDOF(user->quadda, Wloc, &aWquad);CHKERRQ(ierr);
  }

  // need stencil width on locally-computed q
  ierr = DMGetLocalVector(user->quadda,&qloc);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->bloc, &ab);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->m, &am);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da, user->Hinitial, &aHprev);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(user->quadda, qloc, &aqquad);CHKERRQ(ierr);
  // loop over locally-owned elements, including ghosts, to get fluxes at
  // c = 0,1,2,3 points in element;  note start at (xs-1,ys-1)
  for (k = info->ys-1; k < info->ys + info->ym; k++) {
      for (j = info->xs-1; j < info->xs + info->xm; j++) {
          PetscInt  c;
          PetscReal H, Hup;
          Grad      gH, gb;
          for (c=0; c<4; c++) {
              if (user->mtrue) {  // true Mahaffy method
                  // this implementation is at least a factor of two inefficient
                  H = fieldatpt(j,k,locxtrue[c],locytrue[c],aH);
                  gH = gradfatpt(j,k,locxtrue[c],locytrue[c],dx,dy,aH);
                  gb = gradfatpt(j,k,locxtrue[c],locytrue[c],dx,dy,ab);
                  gH = gradav(gH,gradfatpt(j+jnbr[c],k+knbr[c],locxnbr[c],locynbr[c],dx,dy,aH));
                  gb = gradav(gb,gradfatpt(j+jnbr[c],k+knbr[c],locxnbr[c],locynbr[c],dx,dy,ab));
                  Hup = H; // no upwinding allowed in true Mahaffy
              } else { // default M* method
                  H  = fieldatpt(j,k,locx[c],locy[c],aH);
                  gH = gradfatpt(j,k,locx[c],locy[c],dx,dy,aH);
                  gb = gradfatpt(j,k,locx[c],locy[c],dx,dy,ab);
                  if (upwind) {
                      PetscReal lxup = locx[c], lyup = locy[c];
                      if (xdire[c] == PETSC_TRUE)
                          lxup = (gb.x <= 0.0) ? upmin : upmax;
                      else
                          lyup = (gb.y <= 0.0) ? upmin : upmax;
                      Hup = fieldatpt(j,k,lxup,lyup,aH);
                  } else
                      Hup = H;
              }
              if (user->nodiag)
                  aqquad[k][j][c] = getflux(gH,gb,H,Hup,xdire[c],user);
              else
                  aqquad[k][j][c] = getfluxDIAGNOSTIC(gH,gb,H,Hup,xdire[c],user,
                                                      &(aDquad[k][j][c]),&(aWquad[k][j][c]));
          }
      }
  }
  // loop over nodes, not including ghosts, to get residual from quadature over
  // s = 0,1,...,7 points on boundary of control volume (rectangle) around node
  for (k=info->ys; k<info->ys+info->ym; k++) {
      for (j=info->xs; j<info->xs+info->xm; j++) {
          PetscInt s;
          // This is the integral over the control volume boundary using two
          // quadrature points on each side of of the four sides of the
          // rectangular control volume.
          // For M*: two instances of midpoint rule on each side, with two
          //         different values of aq[][] per side
          // For true Mahaffy: ditto, but the two values of aq[][] per side are
          //                   actually the same.
          FF[k][j] = - am[k][j] * dx * dy;
          for (s=0; s<8; s++)
              FF[k][j] += coeff[s] * aqquad[k+ke[s]][j+je[s]][ce[s]];
          if (user->dtres > 0.0)
              FF[k][j] += (aH[k][j] - aHprev[k][j]) * dx * dy / user->dtres;
          if (!user->nodiag) {
              // update diagnostics associated to diffusivity
              PetscReal Dmax = 0.0, Wmagmax = 0.0;
              for (s=0; s<8; s++) {
                  const PetscReal D = aDquad[k+ke[s]][j+je[s]][ce[s]];
                  ds->maxD        = PetscMax(ds->maxD, D);
                  Dmax            = PetscMax(Dmax, D);
                  ds->avD        += D;
                  ds->avDcount   += 1;
                  Wmagmax         = PetscMax(Wmagmax, aWquad[k+ke[s]][j+je[s]][ce[s]]);
              }
              aDnodemax[k][j] = Dmax;
              aWmagnodemax[k][j] = Wmagmax;
          }
      }
  }
  ierr = DMDAVecRestoreArray(user->da, user->m, &am);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, user->Hinitial, &aHprev);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da, user->bloc, &ab);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(user->quadda, qloc, &aqquad);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(user->quadda,&qloc);CHKERRQ(ierr);

  if (!user->nodiag) {
      ierr = DMDAVecRestoreArray(user->da, ds->Dnodemax, &aDnodemax);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArray(user->da, ds->Wmagnodemax, &aWmagnodemax);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArrayDOF(user->quadda, Dloc, &aDquad);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArrayDOF(user->quadda, Wloc, &aWquad);CHKERRQ(ierr);
      ierr = DMRestoreLocalVector(user->quadda,&Dloc);CHKERRQ(ierr);
      ierr = DMRestoreLocalVector(user->quadda,&Wloc);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}


// allows consistent use of "j" and "k" for x,y directions, both in loops and MatSetValuesStencil
typedef struct {
  PetscInt foo,k,j,bar;
} MyStencil;

/* For call-back by SNES using DMDA info.

Evaluates Jacobian matrix on local process patch.

For examples see $PETSC_DIR/src/snes/examples/tutorials/ex5.c or ex9.c.
*/
PetscErrorCode FormJacobianLocal(DMDALocalInfo *info, PetscScalar **aH, Mat jac, Mat jacpre, AppCtx *user) {
  PetscErrorCode  ierr;
  const PetscReal dx = user->dx, dy = dx;
  FLUXINTCOEFFS
  const PetscBool upwind = (user->lambda > 0.0);
  const PetscReal upmin = (1.0 - user->lambda) * 0.5,
                  upmax = (1.0 + user->lambda) * 0.5;
  PetscInt        j, k;
  PetscReal       **ab, ***adQ;
  Vec             dQloc;
  MyStencil       col[33],row;
  PetscReal       val[33];

  PetscFunctionBeginUser;
  if (user->mtrue) {
      SETERRQ(PETSC_COMM_WORLD,1,"ERROR: analytical jacobian not ready in this cases ...\n"); }
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"[inside FormJacobianLocal()]\n"); CHKERRQ(ierr);
  if (user->checkadmissible) {
      ierr = checkadmissible(info,aH); CHKERRQ(ierr); }

  ierr = MatZeroEntries(jac); CHKERRQ(ierr);  // because using ADD_VALUES below

  ierr = DMGetLocalVector(user->sixteenda,&dQloc);CHKERRQ(ierr);

  ierr = DMDAVecGetArray(user->da, user->bloc, &ab);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(user->sixteenda, dQloc, &adQ);CHKERRQ(ierr);
  // loop over locally-owned elements, including ghosts, to get DfluxDl for
  // l=0,1,2,3 at c = 0,1,2,3 points in element;  note start at (xs-1,ys-1)
  for (k = info->ys-1; k < info->ys + info->ym; k++) {
      for (j = info->xs-1; j < info->xs + info->xm; j++) {
          PetscInt  c, l;
          Grad      gH, gb;
          PetscReal H, Hup;
          for (c=0; c<4; c++) {
              PetscReal lxup = locx[c], lyup = locy[c];
              H  = fieldatpt(j,k,locx[c],locy[c],aH);
              gH = gradfatpt(j,k,locx[c],locy[c],dx,dy,aH);
              gb = gradfatpt(j,k,locx[c],locy[c],dx,dy,ab);
              Hup = H;
              if (upwind) {
                  if (xdire[c])
                      lxup = (gb.x <= 0.0) ? upmin : upmax;
                  else
                      lyup = (gb.y <= 0.0) ? upmin : upmax;
                  Hup = fieldatpt(j,k,lxup,lyup,aH);
              }
              for (l=0; l<4; l++) {
                  Grad      dgHdl;
                  PetscReal dHdl, dHupdl;
                  dgHdl  = dgradfatpt(l,j,k,locx[c],locy[c],dx,dy);
                  dHdl   = dfieldatpt(l,j,k,locx[c],locy[c]);
                  dHupdl = (upwind) ? dfieldatpt(l,j,k,lxup,lyup) : dHdl;
                  adQ[k][j][4*c+l] = DfluxDl(gH,gb,dgHdl,H,dHdl,Hup,dHupdl,xdire[c],user);
              }
          }
      }
  }
  // loop over nodes, not including ghosts, to get derivative of residual with respect to nodal value
  for (k=info->ys; k<info->ys+info->ym; k++) {
      row.k = k;
      for (j=info->xs; j<info->xs+info->xm; j++) {
          row.j = j;
          PetscInt s, u, v, l;
          for (s=0; s<8; s++) {
              u = j + je[s];
              v = k + ke[s];
              for (l=0; l<4; l++) {
                  const PetscInt djfroml[4] = { 0,  1,  1,  0},
                                 dkfroml[4] = { 0,  0,  1,  1};
                  col[4*s+l].j = u + djfroml[l];
                  col[4*s+l].k = v + dkfroml[l];
                  val[4*s+l] = coeff[s] * adQ[v][u][4*ce[s]+l];
              }
          }
          if (user->dtjac > 0.0) {
              // add another stencil for diagonal
              col[32].j = j;
              col[32].k = k;
              val[32]   = dx * dy / user->dtjac;
              ierr = MatSetValuesStencil(jac,1,(MatStencil*)&row,33,(MatStencil*)col,val,ADD_VALUES);CHKERRQ(ierr);
          } else {
              ierr = MatSetValuesStencil(jac,1,(MatStencil*)&row,32,(MatStencil*)col,val,ADD_VALUES);CHKERRQ(ierr);
          }
      }
  }
  ierr = DMDAVecRestoreArray(user->da, user->bloc, &ab);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(user->sixteenda, dQloc, &adQ);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(user->sixteenda,&dQloc);CHKERRQ(ierr);

  // Assemble matrix, using the 2-step process:
  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (jacpre != jac) {
    ierr = MatAssemblyBegin(jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode PetscIgnoreZEROPIVOTErrorHandler(MPI_Comm comm,int line,const char *fun,
                   const char *file,PetscErrorCode n,PetscErrorType p,const char *mess,void *ctx) {
   if ((n == PETSC_ERR_MAT_LU_ZRPVT) || (p == PETSC_ERR_MAT_LU_ZRPVT)) {
      int rank;
      MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
      //PetscPrintf(PETSC_COMM_SELF,"---- WARNING: intercepted PETSC_ERR_MAT_LU_ZRPVT on rank %d ----\n",rank);
      AppCtx* user = (AppCtx*)ctx;
      user->luzeropvterr = 1;
      return 0;
   } else {
      return PetscTraceBackErrorHandler(comm,line,fun,file,n,p,mess,ctx);
   }
}

// try a SNES solve and get feedback on the result; H is modified
PetscErrorCode SNESAttempt(SNES *s, Vec H, PetscInt m,
                           SNESConvergedReason *reason, AppCtx *user) {
  PetscErrorCode ierr;
  KSP            ksp;
  PetscInt       its, kspits, luzeropvterr;
  const PetscInt lureason = -99;
  const char     lureasons[30] = "DIVERGED_LU_ZERO_PIVOT";
  char           reasonstr[30];

  PetscFunctionBeginUser;
  user->luzeropvterr = 0;
  PetscPushErrorHandler(PetscIgnoreZEROPIVOTErrorHandler,user);
  SNESSolve(*s, NULL, H);
  PetscPopErrorHandler();
  ierr = MPI_Allreduce(&user->luzeropvterr,&luzeropvterr,1,MPI_INT,MPI_MAX,
                       PETSC_COMM_WORLD); CHKERRQ(ierr);
  if (luzeropvterr > 0)
      *reason = lureason;
  else
      ierr = SNESGetConvergedReason(*s,reason);CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(*s,&its);CHKERRQ(ierr);
  ierr = SNESGetKSP(*s,&ksp); CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp,&kspits); CHKERRQ(ierr);
  strcpy(reasonstr,(*reason==lureason) ? lureasons : SNESConvergedReasons[*reason]);
  myPrintf(user,"%3d. %s   with   ",m,reasonstr);
  myPrintf(user,"eps=%.2e ... %3d KSP (last) iters and %3d Newton iters\n",
           user->eps,kspits,its);
  if (user->dtres > 0.0)
      myPrintf(user,"       (on equations for backward Euler time step of %.3f a)\n",
               user->dtres / user->secpera);
  PetscFunctionReturn(0);
}

