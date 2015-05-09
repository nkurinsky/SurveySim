#!/usr/bin/env python

import os
import sys
import datetime
import time
sys.path.append("../Python/")
from ModelFile import *

mod=ModelFile()
mod.load('test.fits')
mod.info()
mod.update()
