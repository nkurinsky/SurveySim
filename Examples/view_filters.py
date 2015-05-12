#!/usr/bin/env python

#Python prelims
import os
import sys

sys.path.append("../Python/")
from filters import *

instruments=["JWST","Herschel","WISE"]

print(filterDir())

for instrument in instruments:
    ids,names=getFilterIDs(instrument)
    print instrument,names

print(getFilterID("SPIRE_250"))
print(getFilterName(98))
