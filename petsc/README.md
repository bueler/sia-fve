petsc/
======

The code `mahaffy.c` tests classical Mahaffy and improved (M*) finite volume
element methods for the shallow ice approximation (SIA) on a structured Q1 grid.
We use PETSc DMDA and SNESVI to solve the steady SIA free boundary problem.

This code is documented by manuscript `../paper/siafve.tex`.

Build the executable, and run on minimal verification case, by

    $ make mahaffy
    $ ./mahaffy

See comments at start of `mahaffy.c` for dome and bedstep verification runs.
See also `domeconv.sh` and `bedstepconv.sh`.

Some quick verification runs are used as regression tests in

    $ make test

In addition, `mahaffy.c` can read real data and model real ice sheets.  See
the `grn/` and `grn/mcb/` subdirectories for Greenland ice sheet examples.


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

A helpful view of the solution process comes from adding one of

    -snes_monitor_solution -snes_monitor_residual -snes_monitor_solution_update

to the options.  For example, with `-snes_monitor_solution` you see that the
solution `H` starts out too diffuse because the continuation method
over-regularizes at the beginning.

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

