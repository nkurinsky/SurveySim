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

for fname in sys.argv[1:]:
    output=OutputFile(fname)
    model=fname.split('_output')[0]
    field=fname.split('_')[1]
    if(field == "COSMOS"):
        best=0.14
    elif(field == "Lockman-SWIRE"):
        best=0.16
    else:
        raise ValueError("Unknown Field "+field)
    print "{0:<35}{1:.3}".format(model,numpy.exp(best-output.simInfo.fitstat))
