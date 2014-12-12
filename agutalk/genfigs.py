#!/usr/bin/env python
#
# (C) 2014 Ed Bueler
#
# Generate .pdf figures from ASCII VTK written by layer.c:
#   $ (cd ../petsc/ && make layer)
#   $ mkdir foo/
#   $ ../petsc/layer -lay_genfigs foo/
#   $ python genfigs foo/

import numpy
import argparse
import sys

commandline = " ".join(sys.argv[:])

parser = argparse.ArgumentParser(description='FIXME')
# positional filenames
parser.add_argument('prefix', metavar='PREFIX',
                    help='root of file names f-0-0.txt and u-$N-$T.txt',
                    default='')

# process options: simpler names, and type conversions
args = parser.parse_args()
fname = args.prefix + 'f-0-0.txt'
# FIXME
n = 1
t = 0.01
uname = args.prefix + 'u-%d-%g.txt' % (n,t)

# always read .poly file
try:
    ffile = open(fname, 'r')
except IOError:
    print 'cannot open file ', fname
try:
    ufile = open(uname, 'r')
except IOError:
    print 'cannot open file ', uname

headersread = 0
count = 0
for line in ffile:
    if line: # only act if line content remains
        # read headers
        if headersread == 0:
            fheader = line
            #DEBUG print '  header says:  ' + fheader,
            str1 = fheader.split()[0]
            N = (int)(fheader.split()[1])
            print 'reading N = %d values from %s ...' % (N,fname)
            f = numpy.zeros(N)
            headersread += 1
            continue
        elif (headersread == 1) or (headersread == 2):
            fheader = line
            #DEBUG print '  header says:  ' + fheader,
            headersread += 1
            continue
        elif (headersread >= 3) and (count >= N):
            break # nothing more to read
        # read content
        fline = line.split()
        f[count] = float(fline[0])
        #DEBUG print 'f[%3d] = %g' % (count,f[count])
        count += 1


#tikz.write('\\end{tikzpicture}\n')
#tikz.close()
#DEBUG print 'done'

ffile.close()
ufile.close()

