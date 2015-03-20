sia-fve
=======

Project on a stable finite volume element (FVE) method for the steady state
shallow ice approximation (SIA) free boundary problem.

In the paper we construct a modification of the Mahaffy (1976) FD scheme, a new
FVE scheme.  The [PETSc]() code uses a 2D structured rectangular grid, Q1
finite elements, and SNESVI to solve the problem.

This project lived inside my [layer-conserve](https://github.com/bueler/layer-conserve)
repo for much of its life.

