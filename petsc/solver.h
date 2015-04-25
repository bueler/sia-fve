/* (C) 2015 Ed Bueler */

#ifndef SOLVER_H_
#define SOLVER_H_

#include <petscsnes.h>


// indexing of the 8 quadrature points along the boundary of the control volume in M*
// point s=0,...,7 is in element (j,k) = (j+je[s],k+ke[s])
static const PetscInt  je[8] = {0,  0, -1, -1, -1, -1,  0,  0},
                       ke[8] = {0,  0,  0,  0, -1, -1, -1, -1},
                       ce[8] = {0,  3,  1,  0,  2,  1,  3,  2};

// coefficients of quadrature evaluations along the boundary of the control volume in M*
#define FLUXINTCOEFFS \
const PetscReal coeff[8] = {dy/2, dx/2, dx/2, -dy/2, -dy/2, -dx/2, -dx/2, dy/2};

// direction of flux at 4 points in each element
static const PetscBool xdire[4] = {PETSC_TRUE, PETSC_FALSE, PETSC_TRUE, PETSC_FALSE};

// local (element-wise) coords of quadrature points for M*
static const PetscReal locx[4] = {  0.5, 0.75,  0.5, 0.25},
                       locy[4] = { 0.25,  0.5, 0.75,  0.5};

// quadrature points for classical (true) Mahaffy; these are not used in Jacobian
static const PetscReal locxtrue[4] = { 0.5, 1.0, 0.5, 0.0},
                       locytrue[4] = { 0.0, 0.5, 1.0, 0.5},
                       locxnbr[4]  = { 0.5, 0.0, 0.5, 1.0},
                       locynbr[4]  = { 1.0, 0.5, 0.0, 0.5};
static const PetscInt  jnbr[4] = { 0,  1,  0, -1},
                       knbr[4] = {-1,  0,  1,  0};

PetscErrorCode FormFunctionLocal(DMDALocalInfo *info, PetscScalar **aH, PetscScalar **FF, AppCtx *user);

PetscErrorCode FormJacobianLocal(DMDALocalInfo *info, PetscScalar **aH, Mat jac, Mat jacpre, AppCtx *user);

PetscErrorCode SNESAttempt(SNES *s, Vec H, PetscBool again, PetscInt m,
                           SNESConvergedReason *reason, AppCtx *user);


#endif // SOLVER_H_
