petsc/
======

The code `mahaffy.c` tests classical Mahaffy and improved (M*) finite volume
element methods for the shallow ice approximation (SIA) on a structured Q1 grid.
We use PETSc DMDA and SNESVI to solve the steady-state SIA free boundary problem.

This code is documented by the paper in `../paper/`, which was submitted to
the *Journal of Glaciology*.

Quickstart
----------

To build the executable and run a minimal verification case do:

    $ make
    $ ./mahaffy

Some quick verification runs are also used as regression tests:

    $ make test

In addition, `mahaffy.c` can read real data and model real ice sheets.  Also
it can do backward-Euler time-stepping.  See the `grn/`, `grn/robust/`, and
`grn/mcb/` subdirectories for Greenland ice sheet examples.

Usage help and major run modes
------------------------------

Usage help (i.e. options) comes from the usual PETSc `-help` mechanism.  There
are two collections of options, with prefixes `mah_` (for "mahaffy") and `cs_`
(for "continuation scheme"), respectively:

    $ ./mahaffy -help |grep mah_
    $ ./mahaffy -help |grep cs_

Note there are basically these three problem cases.  The first two are
verification cases documented in the paper and in `base/exactsia.{h|c}`.

  1. "Dome" exact solution (the default problem, so `-mah_dome` is optional):

        $ ./mahaffy -mah_dome

  2. "Bedrock step" exact solution:

        $ ./mahaffy -mah_bedstep

  3. Real data from a PETSc binary file.  See `grn/README.md` and
  `grn/mcb/README.md` for two Greenland examples using different bedrock
  topography:

        $ ./mahaffy -mah_read foo.dat

See also `domeconv.sh` and `bedstepconv.sh` for the verification cases.


Link the petsc python scripts
-----------------------------

Both using real data (i.e. using `nc2petsc.py`) and building figures
(`figsmahaffy.py`) requires links to petsc python scripts.  So do this before
proceeding further:

    $ ln -s ~/petsc/bin/PetscBinaryIO.py
    $ ln -s ~/petsc/bin/petsc_conf.py


Figures from dump
-----------------

To make figures, first do this, which writes `unnamed.dat` into `test/`:

    $ mkdir test/
    $ ./mahaffy -mah_dump test/
    $ cd test/

Then run the script which reads the `.dat` file and generates `.png` and
`.pdf` figures:

    $ ../figsmahaffy.py

See `figsmahaffy.py -h` output for additional options.


Initialize with previous result
-------------------------------

The PETSc binary file written with option `-mah_dump` includes the
thickness, so we can read it and restart.  Continuing the previous:

    $ cd ..
    $ ./mahaffy -mah_read test/unnamed.dat -mah_readinitial -cs_start 10


On viewing the solution process
-------------------------------

A helpful graphical view of the solution process comes from adding one of

    -snes_monitor_solution -snes_monitor_residual -snes_monitor_solution_update

to the options.  For example, with `-snes_monitor_solution` you see that the
solution `H` starts out too diffuse because the continuation method
over-regularizes at the beginning.  These options all use X windows, and
often `-draw_pause s`, with `s` in seconds, will be needed.

However, recall we are solving the variational inequality version of equations
`F(H)=0`, with constraint `H>=0`.  With `-snes_monitor_residual` you see that
in the ice-covered area the residual goes to zero as each Newton (SNES)
iteration completes, but that in the _ice-free_ area (i.e. where `H=0`) the
residual `F(H)` is large and positive even when the SNES converges.  This
reflects the complementarity interpretation of the variational inequality, i.e.

    H >= 0  and  F(H) >= 0  and  H F(H) = 0.

Thus the new `-snes_vi_monitor_residual` option in PETSc is helpful because it
shows only the negative values of the residual `F(H)`.

Note `figsmahaffy.py` also shows only the negative values of the residual.

