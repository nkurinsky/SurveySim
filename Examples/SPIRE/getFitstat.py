#!/usr/bin/env python

import os
import sys
import datetime
import time
from OutputFile import *
from pylab import *
import numpy

rc('text', usetex=True)

if (len(sys.argv) < 2):
    print "Please provide model file name"
    quit()

print "{0:<35} {1} {2}".format("Model","ChiSquare","Likelihood")
for fname in sys.argv[1:]:
    output=OutputFile(fname)
    model=fname.split('_output')[0]
    field=fname.split('_')[1]
    if(field == "COSMOS"):
        best=0.1422
    elif(field == "Lockman-SWIRE"):
        best=0.1678
    else:
        raise ValueError("Unknown Field "+field)
    print "{0:<35} {1:.4} {2:.3}".format(model,output.simInfo.fitstat,numpy.exp(best-output.simInfo.fitstat))
