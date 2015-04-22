/* (C) 2015 Ed Bueler */

#ifndef CONTINUATIONSCHEME_H_
#define CONTINUATIONSCHEME_H_

#include <petscvec.h>

#define CSMAX 13
typedef struct {
  PetscInt  max;   // total number of levels:  eps_0, eps_1, ..., eps_{max-1}
  PetscInt  start, // first level to actually use:  eps_start, ..., eps_{end-1}
            end,   // one beyond last level to actually use:  eps_start, ..., eps_{end-1}
            goodm; // last level where convergence happened (not counting recovery)
  PetscReal n0,    // initial value of Glen exponent
            D0,    // initial (and representative?) value of diffusivity
            sched[CSMAX];
} ContinuationScheme;


PetscErrorCode InitializeCS(ContinuationScheme *cs);

PetscErrorCode OptionsCS(ContinuationScheme *cs);

PetscReal epsCS(PetscInt m, const ContinuationScheme *cs);

/* n(1)=n0 and n(0)=n. */
PetscReal nCS(PetscReal n, PetscReal eps, const ContinuationScheme *cs);

/* D(eps)=(1-eps) delta H^{n+2} + eps D_0   so   D(1)=D_0 and D(0)=delta H^{n+2}. */
PetscReal DCS(PetscReal delta, PetscReal H, PetscReal n, PetscReal eps, const ContinuationScheme *cs);

#endif // CONTINUATIONSCHEME_H_
