sia-fve
=======

Project on a stable finite volume element (FVE) method for the steady state
shallow ice approximation (SIA) free boundary problem.

In the `paper/` we construct an improvement on the Mahaffy (1976) FD scheme, a
new FVE scheme.  It uses a 2D structured rectangular grid, Q1 finite elements,
better (over Mahaffy) quadrature in the flux integral, and minimal upwinding.

The C code in `petsc/` uses [PETSc](http://www.mcs.anl.gov/petsc/).
Specifically, a 2d DMDA manages the structured grid SNESVI to solve the problem
which has the non-negative thickness constraint built-in.  Examples include
verification cases and a high-resolution Greenland set-up.

This project lived inside my [layer-conserve](https://github.com/bueler/layer-conserve)
repo for much of its life.

