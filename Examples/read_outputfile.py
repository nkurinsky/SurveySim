#!/usr/bin/env python

import os
import sys
import datetime
import time
sys.path.append("../Python/")
from OutputFile import *

file="spire_output.fits"
output=OutputFile(file)
#output.info()
#output.show()
if(output.fit()):
    output.MCMC.showFit()
    #output.MCMC.showChains()
print ""
