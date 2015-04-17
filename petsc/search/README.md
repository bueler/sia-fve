petsc/search/
=============

We can use the output of `mahaffy` to minimize the "error", i.e. the mis-match
between observed and modeled thickness fields in real-data situations.

Search script
-------------

To test the script use the Rosenbrock objective:

    $ ./searchnA.py --rosen --fprint --show
    $ ./searchnA.py --rosen --nm --disp

For more info:

    $ ./searchnA.py -h

Use mahaffy error as objective
------------------------------

An example in this directory is to search on a grid of values of n (Glen
exponent) and A (ice softness).  First go to directory `../grn/` and preprocess
data, or do `quickstart.sh`, as expected there.  Back in this directory do the
following, which starts with decimating the data to 20km resolution:

    $ cp ../grn/grn.nc .
    $ ncks -d x1,,,4 -d y1,,,4 grn.nc grn20km.nc
    $ ../nc2petsc.py grn20km.nc grn20km.dat

Now run a grid search on 6 processes:

    $ ./searchnA.py --read grn20km.dat --fprint --mpi 6 --show

To use the Nelder-Mead simplex optimization method instead of a grid search, do:

    $ ./searchnA.py --read grn20km.dat --fprint --mpi 6 --disp --nm

To do the above but speed it up with restarting from the last state:

    $ rm -rf foo/
    $ ./searchnA.py --read grn20km.dat --fprint --mpi 6 --disp --nm --restart

This achieves a speed-up by a factor of five, which is significant.


