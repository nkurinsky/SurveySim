#!/usr/bin/env python

#Python prelims
import os
import sys

#codedir=os.getcwd()+'/trunk/';
pydir=os.getcwd()+'/Python/';
seddir=os.getcwd() #codedir+'templates/'
#ensure this is in the path
sys.path.append(pydir)

import time
import pyfits as fits
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import math
from functools import partial
from scipy import interpolate

c = 3.00e08

#hdulist=fits.open(sedfile)
hdu = fits.PrimaryHDU()
hdulist = fits.HDUList([hdu])

lir_list=[9.75,10.0,10.25,10.50,10.75,11.00,11.25,11.50,11.75,12.00,12.25,12.50,12.75,13.00]

sedfile=seddir+'/Kirkpatrick2015_templates.fits'


#============================================================================
#SFR galaxies
#----------------------------------------------------------------------------
with open (seddir+'/Comprehensive_SFG1.txt','r') as f: 
    flines=f.readlines()[4:] #4 lines of header
    fcount=-1
    lam,lnu=[],[]
    for fline in flines:
        tmp1,tmp2,tmp3=fline.split()
        lam.append(float(tmp1)),lnu.append(float(tmp2))

dnu=(np.abs(lam-np.roll(lam,1))/lam)/lam 
dnu=(c*1e6)*dnu
lir=0
for i in range(0,len(lam)-1):
    if(lam[i] >= 5.0):
        lir+=(lnu[i]*dnu[i])/3.86e26

print 'SFG1 template L_TIR', np.log10(lir)
scale=pow(10,lir_list-np.log10(lir))

lam_default=lam
lnu1=np.array(lnu)*scale[0]
lnu2=np.array(lnu)*scale[1]
lnu3=np.array(lnu)*scale[2]
lnu4=np.array(lnu)*scale[3]
lnu5=np.array(lnu)*scale[4]
lnu6=np.array(lnu)*scale[5]
lnu7=np.array(lnu)*scale[6]
lnu8=np.array(lnu)*scale[7]
lnu9=np.array(lnu)*scale[8]

#now the higher luminosity AGN template (defined from 12.1-12.8)
with open (seddir+'/Comprehensive_SFG2.txt','r') as f: 
    flines=f.readlines()[4:]
    fcount=-1
    lam,lnu=[],[]
    for fline in flines:
        tmp1,tmp2,tmp3=fline.split()
        lam.append(float(tmp1)),lnu.append(float(tmp2))

lnu_old=lnu
f = interpolate.interp1d(lam,lnu_old)
#ensure matching boundaries
#print min(lam), max(lam), min(lam_default), max(lam_default)
lnu=f(lam_default)
#these should be the same as above when using the same lambda
#dnu=(np.abs(lam-np.roll(lam,1))/lam)/lam 
#dnu=(c*1e6)*dnu
lir=0
for i in range(0,len(lam_default)-1):
    if(lam_default[i] >= 5.0):
        lir+=(lnu[i]*dnu[i])/3.86e26

print 'SFG2 template L_TIR', np.log10(lir)
scale=pow(10,lir_list-np.log10(lir))

col1=fits.Column(name='lambda1',format='FLOAT',array=lam_default)
col2=fits.Column(name='lnu1',format='FLOAT',array=lnu1)
col3=fits.Column(name='lnu2',format='FLOAT',array=lnu2)
col4=fits.Column(name='lnu3',format='FLOAT',array=lnu3)
col5=fits.Column(name='lnu4',format='FLOAT',array=lnu4)
col6=fits.Column(name='lnu5',format='FLOAT',array=lnu5)
col7=fits.Column(name='lnu6',format='FLOAT',array=lnu6)
col8=fits.Column(name='lnu7',format='FLOAT',array=lnu7)
col9=fits.Column(name='lnu8',format='FLOAT',array=lnu8)
col10=fits.Column(name='lnu9',format='FLOAT',array=lnu9)

col11=fits.Column(name='lnu10',format='FLOAT',array=lnu*scale[9])
col12=fits.Column(name='lnu11',format='FLOAT',array=lnu*scale[10])
col13=fits.Column(name='lnu12',format='FLOAT',array=lnu*scale[11])
col14=fits.Column(name='lnu13',format='FLOAT',array=lnu*scale[12])
col15=fits.Column(name='lnu14',format='FLOAT',array=lnu*scale[13])

cols_SFG_z1=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15])


with open (seddir+'/Comprehensive_SFG3.txt','r') as f: 
    flines=f.readlines()[4:] #4 lines of header
    fcount=-1
    lam,lnu=[],[]
    for fline in flines:
        tmp1,tmp2,tmp3=fline.split()
        lam.append(float(tmp1)),lnu.append(float(tmp2))

dnu=(np.abs(lam-np.roll(lam,1))/lam)/lam 
dnu=(c*1e6)*dnu
lir=0
for i in range(0,len(lam)-1):
    if(lam[i] >= 5.0):
        lir+=(lnu[i]*dnu[i])/3.86e26

print 'SFG3 template L_TIR', np.log10(lir)
scale=pow(10,lir_list-np.log10(lir))

lam_default=lam
lnu1=np.array(lnu)*scale[0]
lnu2=np.array(lnu)*scale[1]
lnu3=np.array(lnu)*scale[2]
lnu4=np.array(lnu)*scale[3]
lnu5=np.array(lnu)*scale[4]
lnu6=np.array(lnu)*scale[5]
lnu7=np.array(lnu)*scale[6]
lnu8=np.array(lnu)*scale[7]
lnu9=np.array(lnu)*scale[8]
lnu10=np.array(lnu)*scale[9]
lnu11=np.array(lnu)*scale[10]
lnu12=np.array(lnu)*scale[11]
lnu13=np.array(lnu)*scale[12]
lnu14=np.array(lnu)*scale[13]

col1=fits.Column(name='lambda1',format='FLOAT',array=lam_default)
col2=fits.Column(name='lnu1',format='FLOAT',array=lnu1)
col3=fits.Column(name='lnu2',format='FLOAT',array=lnu2)
col4=fits.Column(name='lnu3',format='FLOAT',array=lnu3)
col5=fits.Column(name='lnu4',format='FLOAT',array=lnu4)
col6=fits.Column(name='lnu5',format='FLOAT',array=lnu5)
col7=fits.Column(name='lnu6',format='FLOAT',array=lnu6)
col8=fits.Column(name='lnu7',format='FLOAT',array=lnu7)
col9=fits.Column(name='lnu8',format='FLOAT',array=lnu8)
col10=fits.Column(name='lnu9',format='FLOAT',array=lnu9)
col11=fits.Column(name='lnu10',format='FLOAT',array=lnu10)
col12=fits.Column(name='lnu11',format='FLOAT',array=lnu11)
col13=fits.Column(name='lnu12',format='FLOAT',array=lnu12)
col14=fits.Column(name='lnu13',format='FLOAT',array=lnu13)
col15=fits.Column(name='lnu14',format='FLOAT',array=lnu14)

cols_SFG_z2=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15])

#============================================================================
#AGN
#----------------------------------------------------------------------------

with open (seddir+'/Comprehensive_AGN1.txt','r') as f: 
    flines=f.readlines()[13:] #4 lines of header, and 4 lines as the lower-lambda values not well defined for the higher-L templates
    fcount=-1
    lam,lnu=[],[]
    for fline in flines:
        tmp1,tmp2,tmp3=fline.split()
        lam.append(float(tmp1)),lnu.append(float(tmp2))

dnu=(np.abs(lam-np.roll(lam,1))/lam)/lam 
dnu=(c*1e6)*dnu
lir=0
for i in range(0,len(lam)-1):
    if(lam[i] >= 5.0):
        lir+=(lnu[i]*dnu[i])/3.86e26

print 'AGN1 template L_TIR', np.log10(lir)
scale=pow(10,lir_list-np.log10(lir))

#the templates are not on the same scale so will have to rebin the AGN2 to match the AGN1 lambda
lam_default=lam
lnu1=np.array(lnu)*scale[0]
lnu2=np.array(lnu)*scale[1]
lnu3=np.array(lnu)*scale[2]
lnu4=np.array(lnu)*scale[3]
lnu5=np.array(lnu)*scale[4]
lnu6=np.array(lnu)*scale[5]
lnu7=np.array(lnu)*scale[6]
lnu8=np.array(lnu)*scale[7]
lnu9=np.array(lnu)*scale[8]
lnu10=np.array(lnu)*scale[9]


#now the higher luminosity AGN template (defined from 12.1-12.8)
with open (seddir+'/Comprehensive_AGN2.txt','r') as f: 
    flines=f.readlines()[4:]
    fcount=-1
    lam,lnu=[],[]
    for fline in flines:
        tmp1,tmp2,tmp3=fline.split()
        lam.append(float(tmp1)),lnu.append(float(tmp2))

lnu_old=lnu
f = interpolate.interp1d(lam,lnu_old)
#ensure matching boundaries
#print min(lam), max(lam), min(lam_default), max(lam_default)
lnu=f(lam_default)
#these should be the same as above when using the same lambda
#dnu=(np.abs(lam-np.roll(lam,1))/lam)/lam 
#dnu=(c*1e6)*dnu
lir=0
for i in range(0,len(lam_default)-1):
    if(lam_default[i] >= 5.0):
        lir+=(lnu[i]*dnu[i])/3.86e26

print 'AGN2 template L_TIR', np.log10(lir)
scale=pow(10,lir_list-np.log10(lir))

col1=fits.Column(name='lambda1',format='FLOAT',array=lam_default)
col2=fits.Column(name='lnu1',format='FLOAT',array=lnu1)
col3=fits.Column(name='lnu2',format='FLOAT',array=lnu2)
col4=fits.Column(name='lnu3',format='FLOAT',array=lnu3)
col5=fits.Column(name='lnu4',format='FLOAT',array=lnu4)
col6=fits.Column(name='lnu5',format='FLOAT',array=lnu5)
col7=fits.Column(name='lnu6',format='FLOAT',array=lnu6)
col8=fits.Column(name='lnu7',format='FLOAT',array=lnu7)
col9=fits.Column(name='lnu8',format='FLOAT',array=lnu8)
col10=fits.Column(name='lnu9',format='FLOAT',array=lnu9)
col11=fits.Column(name='lnu10',format='FLOAT',array=lnu10)

col12=fits.Column(name='lnu11',format='FLOAT',array=lnu*scale[10])
col13=fits.Column(name='lnu12',format='FLOAT',array=lnu*scale[11])
col14=fits.Column(name='lnu13',format='FLOAT',array=lnu*scale[12])
col15=fits.Column(name='lnu14',format='FLOAT',array=lnu*scale[13])

cols_AGN_z1=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15])


#=========================================================================
# the z-1.8 set
#-------------------------------------------------------------------------

with open (seddir+'/Comprehensive_AGN3.txt','r') as f: 
    flines=f.readlines()[13:] #4 lines of header, and 4 lines as the lower-lambda values not well defined for the higher-L templates
    fcount=-1
    lam,lnu=[],[]
    for fline in flines:
        tmp1,tmp2,tmp3=fline.split()
        lam.append(float(tmp1)),lnu.append(float(tmp2))

dnu=(np.abs(lam-np.roll(lam,1))/lam)/lam 
dnu=(c*1e6)*dnu
lir=0
for i in range(0,len(lam)-1):
    if(lam[i] >= 5.0):
        lir+=(lnu[i]*dnu[i])/3.86e26

print 'AGN3 template L_TIR', np.log10(lir)
scale=pow(10,lir_list-np.log10(lir))

#the templates are not on the same scale so will have to rebin the AGN2 to match the AGN1 lambda
lam_default=lam
lnu1=np.array(lnu)*scale[0]
lnu2=np.array(lnu)*scale[1]
lnu3=np.array(lnu)*scale[2]
lnu4=np.array(lnu)*scale[3]
lnu5=np.array(lnu)*scale[4]
lnu6=np.array(lnu)*scale[5]
lnu7=np.array(lnu)*scale[6]
lnu8=np.array(lnu)*scale[7]
lnu9=np.array(lnu)*scale[8]
lnu10=np.array(lnu)*scale[9]
lnu11=np.array(lnu)*scale[10]
lnu12=np.array(lnu)*scale[11]


#now the higher luminosity AGN template (defined from 12.1-12.8)
with open (seddir+'/Comprehensive_AGN4.txt','r') as f: 
    flines=f.readlines()[4:]
    fcount=-1
    lam,lnu=[],[]
    for fline in flines:
        tmp1,tmp2,tmp3=fline.split()
        lam.append(float(tmp1)),lnu.append(float(tmp2))

lnu_old=lnu
f = interpolate.interp1d(lam,lnu_old)
#ensure matching boundaries
#print min(lam), max(lam), min(lam_default), max(lam_default)
lnu=f(lam_default)
#these should be the same as above when using the same lambda
#dnu=(np.abs(lam-np.roll(lam,1))/lam)/lam 
#dnu=(c*1e6)*dnu
lir=0
for i in range(0,len(lam_default)-1):
    if(lam_default[i] >= 5.0):
        lir+=(lnu[i]*dnu[i])/3.86e26

print 'AGN4 template L_TIR', np.log10(lir)
scale=pow(10,lir_list-np.log10(lir))

col1=fits.Column(name='lambda1',format='FLOAT',array=lam_default)
col2=fits.Column(name='lnu1',format='FLOAT',array=lnu1)
col3=fits.Column(name='lnu2',format='FLOAT',array=lnu2)
col4=fits.Column(name='lnu3',format='FLOAT',array=lnu3)
col5=fits.Column(name='lnu4',format='FLOAT',array=lnu4)
col6=fits.Column(name='lnu5',format='FLOAT',array=lnu5)
col7=fits.Column(name='lnu6',format='FLOAT',array=lnu6)
col8=fits.Column(name='lnu7',format='FLOAT',array=lnu7)
col9=fits.Column(name='lnu8',format='FLOAT',array=lnu8)
col10=fits.Column(name='lnu9',format='FLOAT',array=lnu9)
col11=fits.Column(name='lnu10',format='FLOAT',array=lnu10)

col12=fits.Column(name='lnu11',format='FLOAT',array=lnu11)
col13=fits.Column(name='lnu12',format='FLOAT',array=lnu12)
col14=fits.Column(name='lnu13',format='FLOAT',array=lnu*scale[12])
col15=fits.Column(name='lnu14',format='FLOAT',array=lnu*scale[13])

cols_AGN_z2=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15])



#==========================================================================
# structure the different tables
#---------------------------------------------------------------------------
tbhdu11=fits.new_table(cols_SFG_z1)
tbhdu12=fits.new_table(cols_SFG_z2)

tbhdu21=fits.new_table(cols_AGN_z1)
tbhdu22=fits.new_table(cols_AGN_z2)

tbhdu11.header.update('EXTNAME','SFG_z1',
                        'name of this binary table extension')
tbhdu12.header.update('EXTNAME','SFG_z2',
                        'name of this binary table extension')

tbhdu21.header.update('EXTNAME','AGN_z1',
                        'name of this binary table extension')
tbhdu22.header.update('EXTNAME','AGN_z2',
                        'name of this binary table extension')


tbhdr=tbhdu11.header
tbhdr=tbhdu12.header
tbhdr=tbhdu21.header
tbhdr=tbhdu22.header

hdr=hdulist[0].header
hdulist.close
thdulist=fits.HDUList([hdulist[0],tbhdu11,tbhdu12,tbhdu21,tbhdu22])
if(os.path.isfile(sedfile)):
    print 'Replacing existing template file...'
    os.remove(sedfile)
    thdulist.writeto(sedfile)
else:
    print 'Writing a new template file...'
    thdulist.writeto(sedfile)
