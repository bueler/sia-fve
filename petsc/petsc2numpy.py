#!/usr/bin/env python
# (C) 2015 Ed Bueler

import sys
import numpy as np

try:
    import PetscBinaryIO as pbio
except:
    print "'import PetscBinaryIO' failed"
    print "need link to petsc/bin/petsc-pythonscripts/PetscBinaryIO.py?"
    sys.exit(2)

try:
    import petsc_conf
except:
    print "'import petsc_conf.py' failed"
    print "need link to petsc/bin/petsc-pythonscripts/petsc_conf.py?"
    sys.exit(2)

io = pbio.PetscBinaryIO()

# read Vecs from a specially-formated petsc binary file
# returns a list of numpy arrays
# generally call twice, first to get dimensions, then to get variables
def readvecs(fname,skip=0,num=7,shape=(0,0),failonmissing=True,metaname='variable',names=None):
    try:
        fh = open(fname)
    except:
        if failonmissing:
            print "unable to open '%s' ... ending ..." % fname
            sys.exit(3)
        else:
            return None
    vecs = []
    if names:
        names.reverse()
    for j in range(num):
        try:
            objecttype = io.readObjectType(fh)
        except pbio.DoneWithFile:
            if failonmissing:
                print "no next vec in '%s' ... ending ..." % fname
                sys.exit(2)
            else:
                return vecs
        if objecttype == 'Vec':
            v = io.readVec(fh)
            if (j >= skip):
                if (sum(shape) > 0):
                    v = np.reshape(v,shape)
                vecs.append(v)
                if names:
                    print "  read %s '%s' with shape %s" % (metaname,names.pop(),str(np.shape(v)))
                else:
                    print "  read %s with shape %s" % (metaname,str(np.shape(v)))
        else:
            if failonmissing:
                print "unexpected objecttype '%s' ... ending ..." % objecttype
                sys.exit(4)
            else:
                return None
    fh.close()
    return vecs

