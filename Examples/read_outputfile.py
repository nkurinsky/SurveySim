#!/usr/bin/env python

import os
import sys
import datetime
import time
sys.path.append("../Python/")
from OutputFile import *

files=["herschel_output.fits","herschel_sim_output.fits"]

for file in files:
    output=OutputFile(file)
    output.info()
    output.show()
    if(output.fit()):
        output.MCMC.showFit()
        output.MCMC.showChains()
    print ""
