#!/usr/bin/env python

import os
import sys
import datetime
import time
sys.path.append("/usr/local/surveysim/python/")
from OutputFile import *

if (len(sys.argv) < 2):
    print "Please provide model file name"
else:
    file=sys.argv[1]
    output=OutputFile(file)
    output.info()
    for k,v in output.images.iteritems():
        plt.clf()
        v.plot(cmap=cm.Greys,interpolation='Spline16')
        plt.xlim(1.0,2.0)
        plt.ylim(-2.5,-1.25)
        plt.xlabel("Log(Flux1)")
        plt.ylabel("Color(F1,F2)")
        plt.title('Example Diagnostic: '+k)
        plt.savefig('diag_'+k)
    plt.clf()
    output.saveImages('diags')#,xrange=(0.5,2.5),yrange=(-0.5,0.8))
