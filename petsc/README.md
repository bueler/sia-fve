petsc/
======

The code `mahaffy.c` tests classical Mahaffy and improved (M*) finite volume
element methods for the shallow ice approximation (SIA) on a structured Q1 grid.
We use PETSc DMDA and SNESVI to solve the steady SIA free boundary problem.

This code is documented by manuscript `../paper/mahaffyfem.tex`.

Build the executable, and run on minimal verification case, by

    $ make mahaffy
    $ ./mahaffy

See comments at start of `mahaffy.c` for dome and bedstep verification runs.

In addition, `mahaffy.c` can read real data and model real ice sheets.  See
the `grn/` subdirectory for a Greenland ice sheet example.

Figures from dump
-----------------

To make figures, do:

    $ ln -s ~/petsc/bin/petsc-pythonscripts/PetscBinaryIO.py
    $ ln -s ~/petsc/bin/petsc-pythonscripts/petsc_conf.py
    $ mkdir test/
    $ ./mahaffy -mah_dump test/
    $ cd test/
    $ ../figsmahaffy.py

See `figsmahaffy.py -h` output for additional options.

On viewing the solution process
-------------------------------

A helpful view of the solution process comes from adding one of

    -snes_monitor_solution -snes_monitor_residual -snes_monitor_solution_update

to the options.  With `-snes_monitor_solution` you see that the solution `H`
starts out too diffuse because the continuation method over-regularizes at the
beginning.

However, recall we are solving the variational inequality version of equations
`F(H)=0`, with constraint `H>=0`.  With `-snes_monitor_residual` you see that
in the ice-covered area the residual goes to zero as each Newton (SNES)
iteration completes, but that in the _ice-free_ area (i.e. where `H=0`) the
residual `F(H)` is large and positive even when the SNES converges.  This
reflects the complementarity interpretation of the variational inequality, i.e.

    H >= 0  and  F(H) >= 0  and  H F(H) = 0.

Thus the new `-snes_vi_monitor_residual` option in PETSc (may not be in
maint/master yet) is helpful in the sense of showing only the negative values of
the residual `F(H)`.

The way `figsmahaffy.py` also shows only the negative values of the residual.

