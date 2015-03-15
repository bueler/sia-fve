petsc/mahaffy/
==============

The code `mahaffy.c` tests Mahaffy and Mahaffy* finite volume element methods,
which solve the shallow ice approximation (SIA) on a structured Q1 grid.  It
uses PETSc DMDA and SNESVI to solve the steady SIA free boundary problem.

This code is documented by manuscript `../mahaffy/mahaffyfem.tex`.  It is also
an instance of the general situation documented in `../paper/lc.tex`.

Build the executable, and run on minimal verification case, by

    $ make mahaffy
    $ ./mahaffy

See comments at start of `mahaffy.c` for additional verification runs.

To make figures for the verification results, do:

    $ ln -s ~/petsc/bin/petsc-pythonscripts/PetscBinaryIO.py
    $ ln -s ~/petsc/bin/petsc-pythonscripts/petsc_conf.py
    $ mkdir test/
    $ ./mahaffy -mah_dump test/
    $ cd test/
    $ ../figsmahaffy.py --profile --map

See `figsmahaffy.py -h` output for additional options.

In addition, `mahaffy.c` can read real data and model real ice sheets.  See
the `grn/` subdirectory for a Greenland ice sheet example.
