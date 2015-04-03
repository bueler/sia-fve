petsc/search/
=============

We can use the output of `mahaffy` to minimize the "error", i.e. the mis-match
between observed and modeled thickness fields in real-data situations.

An example in this directory is to search on a grid of values of n (Glen
exponent) and A (ice softness):

    $ cd ..
    $ make mahaffy
    $ cd grn
    $ ncks -d x1,,,4 -d y1,,,4 grn.nc grn20km.nc
    $ ./grn2petsc.py grn20km.nc grn20km.dat
    $ cd search/  # back to here
    $ ./searchnA.py --read ../grn/grn20km.dat --fprint --mpi 6 --show

To use the Nelder-Mead simplex optimization method instead of a grid search, do:

    $ ./searchnA.py --read ../grn/grn20km.dat --fprint --mpi 6 --nm

To simply test the methods use the Rosenbrock objective:

    $ ./searchnA.py --rosen --fprint --show
    $ ./searchnA.py --rosen --fprint --nm --disp

For more info on the search method:

    $ ./searchnA.py -h
