#!/usr/bin/env python

import os
import sys
import datetime
import time
from SurveySim.OutputFile import *
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
    print model
    output.MCMC.info()
    print ""
