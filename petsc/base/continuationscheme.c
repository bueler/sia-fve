/* (C) 2015 Ed Bueler */

#include "continuationscheme.h"

PetscErrorCode SetFromOptionsCS(const char *optprefix, ContinuationScheme *cs) {
  PetscErrorCode ierr;
  PetscInt i;
  cs->max   = CSMAX;
  cs->start = 0;
  cs->end   = cs->max;
  for (i = 0; i < cs->max - 1; i++)
      cs->sched[i] = PetscPowReal(0.1,((double)i) / 3.0);
  cs->sched[cs->max-1] = 0.0;
  cs->n0    = 1.0;
  cs->D0    = 10.0;        // m^2 / s
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,optprefix,"options to continuation scheme","");CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-D0", "initial (and representative?) value of diffusivity; in m^2/s",
      NULL,cs->D0,&cs->D0,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt(
      "-end", "one beyond last level to actually use:  eps_start, ..., eps_{end-1}",
      NULL,cs->end,&cs->end,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal(
      "-n0", "initial value of Glen exponent",
      NULL,cs->n0,&cs->n0,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt(
      "-start", "first level to actually use:  eps_start, ..., eps_{end-1}",
      NULL,cs->start,&cs->start,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  if (cs->end > cs->max) {
      SETERRQ(PETSC_COMM_WORLD,1,"ERROR in ContinuationScheme option: -cs_end > cs.max\n");
  }
  if (cs->start >= cs->end) {
      SETERRQ(PETSC_COMM_WORLD,1,"ERROR in ContinuationScheme option: -cs_start >= -cs_end\n");
  }
  if (cs->start < 0) {
      SETERRQ(PETSC_COMM_WORLD,1,"ERROR in ContinuationScheme option: -cs_start < 0\n");
  }
  if (cs->end < 1) {
      SETERRQ(PETSC_COMM_WORLD,1,"ERROR in ContinuationScheme option: -cs_end < 1\n");
  }
  if (cs->D0 <= 0.0) {
      SETERRQ(PETSC_COMM_WORLD,1,"ERROR in ContinuationScheme option: -cs_D0 <= 0.0\n");
  }
  if (cs->n0 < 1.0) {
      SETERRQ(PETSC_COMM_WORLD,1,"ERROR in ContinuationScheme option: -cs_n0 < 1.0\n");
  }
  PetscFunctionReturn(0);
}

PetscInt startCS(const ContinuationScheme *cs) {
  return cs->start;
}

PetscInt endCS(const ContinuationScheme *cs) {
  return cs->end;
}

PetscReal epsCS(PetscInt m, const ContinuationScheme *cs) {
  return cs->sched[m];
}

/* n(1)=n0 and n(0)=n. */
PetscReal nCS(PetscReal n, PetscReal eps, const ContinuationScheme *cs) {
  return (1.0 - eps) * n + eps * cs->n0;
}

/* D(eps)=(1-eps) delta H^{n+2} + eps D_0   so   D(1)=D_0 and D(0)=delta H^{n+2}. */
PetscReal DCS(PetscReal delta, PetscReal H, PetscReal n, PetscReal eps, const ContinuationScheme *cs) {
  return (1.0 - eps) * delta * PetscPowReal(PetscAbsReal(H),n+2.0) + eps * cs->D0;
}

