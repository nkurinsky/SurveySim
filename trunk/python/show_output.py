#!/usr/bin/env python

import os
import sys
import datetime
import time
#sys.path.append("/usr/local/surveysim/python/")
from OutputFile import *

if (len(sys.argv) < 2):
    print "Please provide model file name"
else:
    file=sys.argv[1]
    output=OutputFile(file)
    output.info()
#    output.saveCounts('spire_counts.pdf')
    output.show()
#    if(output.fit()):
#        output.MCMC.showFit()
#        output.MCMC.showChains()
    print ""
