paper/
======

The manuscript shows that the Mahaffy method is an FVE method, and then improves
it by better quadrature and selective upwinding.  The parallel, Newton
implementation solves it as a free boundary problem in steady state.  Thus we
see high accuracy, including in non-flat bed applicaiton, and no stability
restrictions whatsoever.

To build `mahaffyfem.pdf`, the IGS classes are needed, specifically `igs.cls`,
`igs.bst`, and `igsnatbib.sty` fro the v2 version of the IGS Latex class.

Then do

    $ make

