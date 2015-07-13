#!/usr/bin/python

## Debugging filter.py
## Problem: getFilterID() is raising a LookupError.
##          It can't find AzTEC 1.1 is FILTER_LIST.
##          Problem is in read_filters()
##          read_filters() is not reading the newly added line,
##          It says that there are only 125 items in FILTER_LIST,
##          when there should be 126.

import os
import numpy as np
from filters import *

fdir = filterDir()
ids,names = read_filters()

print len(ids)
print len(names)

print ids



## Solved: filterDir() was grabbing old surveysim directory 
