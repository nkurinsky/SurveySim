#!/usr/bin/env python

import os
import sys
import datetime
import time
sys.path.append('/Users/Jed/SurveySim/trunk/python/')
from OutputFile import *

if (len(sys.argv) < 2):
    print 'Please provide model file name'
else:
    file=sys.argv[1]
    output=OutputFile(file)
    output.info()
