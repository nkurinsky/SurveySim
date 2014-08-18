#!/usr/bin/python
#
# Used to plot new filters in this directory to ensure compatibility with program or view filter profile easily
#

import sys
import numpy as np
from pylab import *

def main(argv):
    if(len(argv) < 1):
        print("Usage: plot_filter.py filename")
        exit(1)
        
    print("Filters:")
    print(argv)
    clf()
    for arg in argv:
        w,f = np.loadtxt(arg,unpack=True,usecols=[0,1])
        plot(log(w),f,label=arg[:len(arg)-4])
    xlabel("Log(Wavelength [microns])")
    ylabel("Spectral Response Curve")
    ylim(0,1.2)
    legend(loc='upper center',ncol=len(argv)/4+1,prop={'size':12})
    show()
    exit(0)

if __name__ == "__main__":
    main(sys.argv[1:])

