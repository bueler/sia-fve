/* (C) 2015 Ed Bueler */

#ifndef MAHAFFYCTX_H_   /* Include guard */
#define MAHAFFYCTX_H_

#include <petscdmda.h>
#include <petscvec.h>

typedef struct {
  DM        da, quadda;
  Vec       b,      // the bed elevation
            m,      // the (steady) surface mass balance
            Hexact; // the exact thickness (either in verification or data sense)
  PetscReal dx,     // fixed grid spacing; dx = dy
            Lx,     // domain is [-Lx,Lx] x [-Ly,Ly]
            Ly,
            n,      // Glen exponent for SIA flux term
            g,      // acceleration of gravity
            rho,    // ice density
            secpera,// number of seconds in a year
            A,      // ice softness
            Gamma,  // coefficient for SIA flux term
            eps,    // dimensionless regularization
            maxD,   // local value used for reporting maximum of diffusivity
            D0;     // representative value of diffusivity (in regularization)
  PetscInt  Nx,     // grid has Nx x Ny nodes
            Ny,
            Neps;   // number of levels in regularization/continuation
  PetscBool mtrue,  // use true Mahaffy method instead of Mahaffy* (default)
            read,   // read grid and data from PETSc binary file
            dome,   // use dome exact solution ("Bueler profile")
            bedstep,// use bedrock step exact solution from Jarosch, Schoof, Anslow (2013)
            swapxy, // swap x and y axes in building exact solution
            showdata,// show b and m with X viewer
            dump;   // dump fields into ASCII VTK files
  char      figsprefix[512],
            readname[512];
} AppCtx;

#endif // MAHAFFYCTX_H_
