import pyfits
import numpy
from pylab import *

hdus = pyfits.open("spire_fls_dr2.fits")
ftable = hdus[1].data

f1 = ftable["F250"]
ef1 = ftable["EF250"]

f2 = ftable["F350"]
ef2 = ftable["EF350"]

f3 = ftable["F500"]
ef3 = ftable["EF500"]

gpts1 = numpy.where(f1 > 0)
gpts2 = numpy.where(f2 > 0)
gpts3 = numpy.where(f3 > 0)

print(len(gpts1[0]))
print(len(gpts2[0]))
print(len(gpts3[0]))
