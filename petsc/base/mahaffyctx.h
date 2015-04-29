/* (C) 2015 Ed Bueler */

/* Structs which are used by mahaffy.c. */

#ifndef MAHAFFYCTX_H_
#define MAHAFFYCTX_H_

#include <sys/time.h>
#include <petscdmda.h>
#include <petscvec.h>

#include "continuationscheme.h"

// holds value of gradient at point
typedef struct {
    PetscReal x,y;
} Grad;

typedef struct {
  Vec       Dnodemax,// maximum value of the diffusivity D at the quadrature points for that node
            Wmagnodemax;// maximum value of the magnitude of the pseudo-velocity W at the
                        // quadrature points for that node
  PetscReal avD,    // local average value of diffusivity, in m^2 s^-1
            maxD;   // local maximum value of diffusivity, in m^2 s^-1
  PetscInt  avDcount;// used to get local average value of diffusivity
} DiagnosticScheme;

typedef struct {
  DM        da, quadda, sixteenda;
  Vec       b,      // the bed elevation
            m,      // the (steady) surface mass balance
            Hexact, // the exact or observed thickness (verification or data, resp.)
            Hinitial,// holds initial state (if read from file) and used when doing time step
            bloc;   // copy of bed elevation with ghosts
  PetscReal dx,     // fixed grid spacing; dx = dy
            Lx,     // domain is [-Lx,Lx] x [-Ly,Ly]
            Ly,
            n,      // Glen exponent for SIA flux term
            g,      // acceleration of gravity
            rho,    // ice density
            secpera,// number of seconds in a year
            A,      // ice softness
            Gamma,  // coefficient for SIA flux term
            eps,    // current continuation parameter for n and D
            delta,  // dimensionless regularization for slope in SIA formulas
            lambda, // amount of upwinding; lambda=0 is none and lambda=1 is "full"
            dtres,  // time step for backward Euler used in residual
            dtjac,  //     ...                      used in Jacobian
            dtrecovery,//  ...                      used in recovery
            initmagic;// constant, in years, used to multiply SMB to get initial iterate for thickness
  PetscInt  Nx, Ny, // grid has Nx x Ny nodes
            numsteps,// number of time-steps to take; is set to 1 if steady-state
            recoverycount,// number of steps of recovery taken; zero if recovery did not happen
            luzeropvterr, // error handler sets this if it "intercepts" zero pivot error
            goodm;  // last continuation scheme level where convergence happened (not counting recovery)
  PetscBool mtrue,  // use true Mahaffy method instead of M*
            read,   // read grid and data from special-format PETSc binary file
            readinitial,// use read initial H instead of generating guess in usual way
            dome,   // use dome exact solution ("Bueler profile")
            bedstep,// use bedrock step exact solution from Jarosch, Schoof, Anslow (2013)
            swapxy, // swap x and y axes in building exact solution
            showdata,// show b and m with X viewer
            checkadmissible,// in FormFunctionLocal(), stop if H < 0.0
            divergetryagain,// on SNES diverge, try again with eps *= 1.5
            dump,   // dump fields into individual PETSc binary files
            silent, // run silent
            averr,  // only display average error at end
            maxerr, // only display maximum error at end
            nodiag, // do not use DiagnosticScheme
            history;// write ASCII history file
  char      figsprefix[512],
            readname[512],
            readinitialname[512];
  struct timeval     starttime, endtime;
  ContinuationScheme cs;
  DiagnosticScheme   ds;
} AppCtx;

#endif // MAHAFFYCTX_H_

