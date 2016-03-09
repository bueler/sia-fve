paper/
======

This paper shows that the Mahaffy method is an FVE method, and then improves
it by better quadrature and selective upwinding.  The parallel, Newton
implementation solves it as a free boundary problem in steady state.  Thus we
see high accuracy, including in non-flat bed application, and no stability
restrictions whatsoever.

It was submitted 5/17/2015 to the *Journal of Glaciology* and accepted
10/15/2015.  As of 3/9/2016 it is still "to appear".

To build `siafve.pdf`, the IGS classes are needed, specifically `igs.cls`,
`igs.bst`, `igsnatbib.sty`, and `lineno.sty` from the v3.02 version
of the IGS Latex class.

Then do

    $ make

