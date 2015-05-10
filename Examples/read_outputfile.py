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
    print file
    print output.type()
    print output.fit()
    output.info()
    output.parameters.info()
    if(output.fit()):
        plt.figure(figsize=(12,5))
        output.showImages()
        output.MCMC.showR()
        output.MCMC.showFit()
    print ""
