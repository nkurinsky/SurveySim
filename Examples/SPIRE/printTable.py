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

possKeys=['L0','PHI0', 'ALPHA', 'BETA', 'P', 'Q', 'P2', 'Q2', 'ZBP', 'ZBQ', 'FAGN0', 'T1', 'T2', 'ZBT', 'FCOMP', 'FCOLD']

headStr="{0:<35}\t\t{1}".format("Model","Chi-Sq")
for key in possKeys:
    headStr=headStr+"\t{0}".format(key)
print headStr
for fname in sys.argv[1:]:
    output=OutputFile(fname)
    model=fname.split('_output')[0]
    field=fname.split('_')[1]
    if(field == "COSMOS"):
        best=0.14224
    elif(field == "Lockman-SWIRE"):
        best=0.16784
    else:
        raise ValueError("Unknown Field "+field)

    pnames=output.simInfo.initial.keys()
    outString="{0:<35}\t\t{1:.4}".format(model,output.simInfo.fitstat)
    for key in possKeys:
        if key in pnames:
            if(output.simInfo.initial[key] == 0.0):
                outString=outString+"\t  -  "
            else:
                outString=outString+"\t{0:.3}".format(output.simInfo.initial[key])
        else:
             outString=outString+"\t  -  "

    print outString
