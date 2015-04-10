#!/usr/bin/env python

import pyfits

lims1=9.3
lims2=9.6
lims3=13.5
error1=3.1
error2=3.2
error3=4.5

hdus=pyfits.open("trunk/model/model.fits",mode='update')

hdus[0].header['ERROR1']=error1
hdus[0].header['ERROR2']=error2
hdus[0].header['ERROR3']=error3
hdus[0].header['LIMIT1']=lims1
hdus[0].header['LIMIT2']=lims2
hdus[0].header['LIMIT3']=lims3

hdus.flush()
hdus.close()
