/* (C) 2015 Ed Bueler */

/* Structs which are used by mahaffy.c. */

#ifndef MAHAFFYCTX_H_
#define MAHAFFYCTX_H_

#include <sys/time.h>
#include <petscdmda.h>
#include <petscvec.h>

// holds value of gradient or vector at point
typedef struct {
    PetscReal x,y;
} Grad;

typedef struct {
  DM        da, quadda, sixteenda;
  Vec       b,      // the bed elevation
            m,      // the (steady) surface mass balance
            Hexact, // the exact or observed thickness (verification or data, resp.)
            Hprev,  // only used when doing backward Euler time step as recovery
            bloc,   // copy of bed elevation with ghosts
            Dquad,  // diffusivity D at the quadrature points on an element
            Dnodemax,// maximum value of the diffusivity D at the quadrature points for that node
            Wmagnodemax;// maximum value of the magnitude of the pseudo-velocity W at the
                        // quadrature points for that node
  PetscReal dx,     // fixed grid spacing; dx = dy
            Lx,     // domain is [-Lx,Lx] x [-Ly,Ly]
            Ly,
            n,      // Glen exponent for SIA flux term
            g,      // acceleration of gravity
            rho,    // ice density
            secpera,// number of seconds in a year
            A,      // ice softness
            Gamma,  // coefficient for SIA flux term
            eps,    // current dimensionless regularization for n and D
            delta,  // current dimensionless regularization for slope in SIA formulas
            lambda, // amount of upwinding; lambda=0 is none and lambda=1 is "full"
            dtBE,   // time step for backward Euler, used in recovery
            avD,    // local average value of diffusivity, in m^2 s^-1; used for reporting
            maxD,   // local maximum value of diffusivity, in m^2 s^-1; used for reporting
            n0,     // initial value of Glen exponent in continuation scheme
            D0,     // initial (and representative?) value of diffusivity in continuation scheme
            initmagic,// constant, in years, used to multiply SMB to get initial iterate for thickness
            eps_sched[13];
  PetscInt  Nx,     // grid has Nx x Ny nodes
            Ny,
            Neps,   // number of levels in continuation scheme
            avDcount,// used to get local average value of diffusivity; used for reporting
            luzeropvterr; // error handler sets this if it "intercepts" zero pivot error
  PetscBool mtrue,  // use true Mahaffy method instead of M*
            read,   // read grid and data from special-format PETSc binary file
            dome,   // use dome exact solution ("Bueler profile")
            bedstep,// use bedrock step exact solution from Jarosch, Schoof, Anslow (2013)
            swapxy, // swap x and y axes in building exact solution
            showdata,// show b and m with X viewer
            checkadmissible,// in FormFunctionLocal(), stop if H < 0.0
            divergetryagain,// on SNES diverge, try again with eps *= 1.5
            dorecovery,// on SNES diverge, set this to TRUE for recovery attempt (by whatever mechanism)
            dump,   // dump fields into individual PETSc binary files
            silent, // run silent
            averr,  // only display average error at end
            maxerr, // only display maximum error at end
            nodiag, // do not store, generate, or output diagnostic D and Wmag fields
            history;// write ASCII history file
  char      figsprefix[512],
            readname[512];
  struct timeval starttime,
                 endtime;
} AppCtx;

#endif // MAHAFFYCTX_H_
