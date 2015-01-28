layer-conserve/petsc/
==============

The code `layer.c` describes a one-dimensional moving layer with a non-negative-
constrained thickness, and thus a moving boundary.  PETSc SNESVI is used to
solve the free boundary problem at each time step.  For more information build
`doclayer.pdf` by

    $ make doclayer.pdf

Building `layer.c` requires PETSc 3.5.2 or maint.  Do

    $ make layer

To thoroughly test with exact solution,

    $ ./convtest.sh
