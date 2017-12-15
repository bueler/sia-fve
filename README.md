sia-fve
=======

A stable, implicit finite volume element (FVE) method for the steady state shallow ice approximation (SIA) free-boundary problem.  Published as:

E. Bueler (2016).  _Stable finite volume element schemes for the shallow-ice approximation_, J. Glaciol. 62 (232), 230-242, [doi:10.1017/jog.2015.3](http://dx.doi.org/10.1017/jog.2015.3).

We re-interpret the classical Mahaffy (1976) FD scheme as an FVE scheme.  Then we construct an improved scheme.  Both the classical and improved schemes use a 2D structured rectangular grid, piecewise-bilinear (i.e. Q1 finite element) trial function space, and a flux integral as the weak form.  The improved scheme has better quadrature in the flux integral and a form of first-order upwinding which only acts on the bedrock-gradient part of the flux.  We solve the steady-state problem by a continuation-modified and constrained Newton solver.

Directory `paper/` contains the LaTeX sources and figures.

The C code in `petsc/` implements both schemes using [PETSc](http://www.mcs.anl.gov/petsc/).  Specifically, a 2d DMDA manages the structured grid SNESVI to solve the problem.  The free boundary problem thus has the non-negative thickness constraint built-in.  Examples include verification cases and two high-resolution Greenland set-ups.

This project lived inside my [layer-conserve](https://github.com/bueler/layer-conserve) repo for much of its life.
