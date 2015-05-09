#!/usr/bin/env python

#Python prelims
import os
import sys

sys.path.append("../Python/")
from filters import *

instruments=["JWST","SPIRE","WISE"]

print(filterDir())

for instrument in instruments:
    ids,names=getFilterIDs(instrument)
    print(names)

print(getFilterID("SPIRE_250"))
print(getFilterName(98))
