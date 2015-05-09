#!/usr/bin/env python

import os
import sys
import datetime
import time
sys.path.append("../Python/")
from ModelFile import *

mod=ModelFile()
mod.filters[0].setID("SPIRE_250")
mod.filters[1].setID("SPIRE_350")
mod.filters[2].setID("SPIRE_500")

mod.write('newfile.fits')
mod.info()
