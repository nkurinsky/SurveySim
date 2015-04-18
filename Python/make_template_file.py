#!/usr/bin/env python

#Python prelims
import os
import sys

#codedir=os.getcwd()+'/trunk/';
pydir=os.getcwd();
codedir = pydir.replace(' ', '')[:-6] #remove last 6 characters (i.e. Python)
seddir=codedir+'trunk/templates/'

#ensure this is in the path
#sys.path.append(pydir)

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
#use Rieke+2008 for z=0 and Kirkpatrick+2015 for z=1 and z=2)
#define a common lambda array which goes from 0.5-1000um
#adopt same luminosity bins as in Rieke et al. i.e. 9.75-13.0 in L_TIR(5-1000um)
#----------------------------------------------------------------------------
with open (seddir+'/Rieke2009_templates.txt','r') as f:
    flines=f.readlines()[59:539]
    fcount=-1
    lam_rieke,lnu1,lnu2,lnu3,lnu4,lnu5,lnu6,lnu7,lnu8,lnu9,lnu10,lnu11,lnu12,lnu13,lnu14=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    for fline in flines:
        tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12,tmp13,tmp14=fline.split()
        lam_rieke.append(float(tmp0)),lnu1.append(float(tmp1)),lnu2.append(float(tmp2)),lnu3.append(float(tmp3)),lnu4.append(float(tmp4)),lnu5.append(float(tmp5)),lnu6.append(float(tmp6)),lnu7.append(float(tmp7)),lnu8.append(float(tmp8)),lnu9.append(float(tmp9)),lnu10.append(float(tmp10)),lnu11.append(float(tmp11)),lnu12.append(float(tmp12)),lnu13.append(float(tmp13)),lnu14.append(float(tmp14))

#there are in Jy normalized to D_L=10Mpc
# so to convert to luminosity in W/Hz need to use a conversion factor of
#10^{-26}*4*!pi*(10*10^6*3*10^16)^2.
conv_whz=4*np.pi*9*1e18

lnu1=np.array(lnu1)*conv_whz
lnu2=np.array(lnu2)*conv_whz
lnu3=np.array(lnu3)*conv_whz
lnu4=np.array(lnu4)*conv_whz
lnu5=np.array(lnu5)*conv_whz
lnu6=np.array(lnu6)*conv_whz
lnu7=np.array(lnu7)*conv_whz
lnu8=np.array(lnu8)*conv_whz
lnu9=np.array(lnu9)*conv_whz
lnu10=np.array(lnu10)*conv_whz
lnu11=np.array(lnu11)*conv_whz
lnu12=np.array(lnu12)*conv_whz
lnu13=np.array(lnu13)*conv_whz
lnu14=np.array(lnu14)*conv_whz

#want to extend to 0.5um 
lam_short=[0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.10,4.20,4.30,4.40,4.50,4.60,4.70,5.0,5.1,5.2]
lam_short=np.array(lam_short)
lam_rieke=np.array(lam_rieke)
lam_default=np.concatenate((lam_short,lam_rieke),axis=0)

dnu_default=(np.abs(lam_default-np.roll(lam_default,1))/lam_default)/lam_default 
dnu_default=(c*1e6)*dnu_default

#read in Sc spiral galaxy template from SWIRE template library in order to get the stellar bump and switch to the above at 0.5-4um
with open ('/Users/annie/code/idl/auxdata/swire/Sc_template_norm.sed','r') as f:
    flines=f.readlines()
    fcount=-1
    lam_swire,flam_swire=[],[]
    for fline in flines:
        tmp0,tmp1=fline.split()
        lam_swire.append(float(tmp0)),flam_swire.append(float(tmp1))

lam_swire=np.array(lam_swire)/1e4 #convert to microns
lnu_swire=np.array(flam_swire)*lam_swire**2. #convert to fnu but overall normalization not set yet!
 
f = interpolate.interp1d(lam_swire,lnu_swire)
lnu_tmp=f(lam_short)
#print lam_short[44],lnu_tmp[44],lam_rieke[0],lnu1[0]
tmp=lnu_tmp*(lnu1[0]/lnu_tmp[44])
lnu1_default=np.concatenate((tmp,lnu1),axis=0)
tmp=lnu_tmp*(lnu2[0]/lnu_tmp[44])
lnu2_default=np.concatenate((tmp,lnu2),axis=0)
tmp=lnu_tmp*(lnu3[0]/lnu_tmp[44])
lnu3_default=np.concatenate((tmp,lnu3),axis=0)
tmp=lnu_tmp*(lnu4[0]/lnu_tmp[44])
lnu4_default=np.concatenate((tmp,lnu4),axis=0)
tmp=lnu_tmp*(lnu5[0]/lnu_tmp[44])
lnu5_default=np.concatenate((tmp,lnu5),axis=0)
tmp=lnu_tmp*(lnu6[0]/lnu_tmp[44])
lnu6_default=np.concatenate((tmp,lnu6),axis=0)
tmp=lnu_tmp*(lnu7[0]/lnu_tmp[44])
lnu7_default=np.concatenate((tmp,lnu7),axis=0)
tmp=lnu_tmp*(lnu8[0]/lnu_tmp[44])
lnu8_default=np.concatenate((tmp,lnu8),axis=0)
tmp=lnu_tmp*(lnu9[0]/lnu_tmp[44])
lnu9_default=np.concatenate((tmp,lnu9),axis=0)
tmp=lnu_tmp*(lnu10[0]/lnu_tmp[44])
lnu10_default=np.concatenate((tmp,lnu10),axis=0)
tmp=lnu_tmp*(lnu11[0]/lnu_tmp[44])
lnu11_default=np.concatenate((tmp,lnu11),axis=0)
tmp=lnu_tmp*(lnu12[0]/lnu_tmp[44])
lnu12_default=np.concatenate((tmp,lnu12),axis=0)
tmp=lnu_tmp*(lnu13[0]/lnu_tmp[44])
lnu13_default=np.concatenate((tmp,lnu13),axis=0)
tmp=lnu_tmp*(lnu14[0]/lnu_tmp[44])
lnu14_default=np.concatenate((tmp,lnu14),axis=0)

col1=fits.Column(name='lambda',format='FLOAT',array=lam_default)
col2=fits.Column(name='lnu1',format='FLOAT',array=lnu1_default)
col3=fits.Column(name='lnu2',format='FLOAT',array=lnu2_default)
col4=fits.Column(name='lnu3',format='FLOAT',array=lnu3_default)
col5=fits.Column(name='lnu4',format='FLOAT',array=lnu4_default)
col6=fits.Column(name='lnu5',format='FLOAT',array=lnu5_default)
col7=fits.Column(name='lnu6',format='FLOAT',array=lnu6_default)
col8=fits.Column(name='lnu7',format='FLOAT',array=lnu7_default)
col9=fits.Column(name='lnu8',format='FLOAT',array=lnu8_default)
col10=fits.Column(name='lnu9',format='FLOAT',array=lnu9_default)
col11=fits.Column(name='lnu10',format='FLOAT',array=lnu10_default)
col12=fits.Column(name='lnu11',format='FLOAT',array=lnu11_default)
col13=fits.Column(name='lnu12',format='FLOAT',array=lnu12_default)
col14=fits.Column(name='lnu13',format='FLOAT',array=lnu13_default)
col15=fits.Column(name='lnu14',format='FLOAT',array=lnu14_default)

#plt.plot(lam_default,lnu4_default)
#plt.xscale('log')
#plt.yscale('log')
#plt.show()

cols_SFG_z0=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15])


#read-in the Kirkpatrick et al. SFG templates
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

lam_old=np.array(lam)
lnu_old=np.array(lnu)
if(lam_old[0] > lam_default[0]):
    lam_old=np.concatenate(([0.5],lam_old,[1031]),axis=0)
    lnu_old=np.concatenate(([0],lnu_old,[0]),axis=0)
else:
    lam_old=np.concatenate((lam_old,[1031]),axis=0)
    lnu_old=np.concatenate((lnu_old,[0]),axis=0)

f = interpolate.interp1d(lam_old,lnu_old)
lnu=f(lam_default)

col2=fits.Column(name='lnu1',format='FLOAT',array=lnu*scale[0])
col3=fits.Column(name='lnu2',format='FLOAT',array=lnu*scale[1])
col4=fits.Column(name='lnu3',format='FLOAT',array=lnu*scale[2])
col5=fits.Column(name='lnu4',format='FLOAT',array=lnu*scale[3])
col6=fits.Column(name='lnu5',format='FLOAT',array=lnu*scale[4])
col7=fits.Column(name='lnu6',format='FLOAT',array=lnu*scale[5])
col8=fits.Column(name='lnu7',format='FLOAT',array=lnu*scale[6])
col9=fits.Column(name='lnu8',format='FLOAT',array=lnu*scale[7])
col10=fits.Column(name='lnu9',format='FLOAT',array=lnu*scale[8])

with open (seddir+'/Comprehensive_SFG2.txt','r') as f: 
    flines=f.readlines()[4:]
    fcount=-1
    lam,lnu=[],[]
    for fline in flines:
        tmp1,tmp2,tmp3=fline.split()
        lam.append(float(tmp1)),lnu.append(float(tmp2))

lam_old=lam
lnu_old=lnu
if(lam_old[0] > lam_default[0]):
    lam_old=np.concatenate(([0.5],lam_old,[1031]),axis=0)
    lnu_old=np.concatenate(([0],lnu_old,[0]),axis=0)
else:
    lam_old=np.concatenate((lam_old,[1031]),axis=0)
    lnu_old=np.concatenate((lnu_old,[0]),axis=0)

f = interpolate.interp1d(lam_old,lnu_old)
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

lam_old=np.array(lam)
lnu_old=np.array(lnu)
if(lam_old[0] > lam_default[0]):
    lam_old=np.concatenate(([0.5],lam_old,[1031]),axis=0)
    lnu_old=np.concatenate(([0],lnu_old,[0]),axis=0)
else:
    lam_old=np.concatenate((lam_old,[1031]),axis=0)
    lnu_old=np.concatenate((lnu_old,[0]),axis=0)

f = interpolate.interp1d(lam_old,lnu_old)
lnu=f(lam_default)

col2=fits.Column(name='lnu1',format='FLOAT',array=lnu*scale[0])
col3=fits.Column(name='lnu2',format='FLOAT',array=lnu*scale[1])
col4=fits.Column(name='lnu3',format='FLOAT',array=lnu*scale[2])
col5=fits.Column(name='lnu4',format='FLOAT',array=lnu*scale[3])
col6=fits.Column(name='lnu5',format='FLOAT',array=lnu*scale[4])
col7=fits.Column(name='lnu6',format='FLOAT',array=lnu*scale[5])
col8=fits.Column(name='lnu7',format='FLOAT',array=lnu*scale[6])
col9=fits.Column(name='lnu8',format='FLOAT',array=lnu*scale[7])
col10=fits.Column(name='lnu9',format='FLOAT',array=lnu*scale[8])
col11=fits.Column(name='lnu10',format='FLOAT',array=lnu*scale[9])
col12=fits.Column(name='lnu11',format='FLOAT',array=lnu*scale[10])
col13=fits.Column(name='lnu12',format='FLOAT',array=lnu*scale[11])
col14=fits.Column(name='lnu13',format='FLOAT',array=lnu*scale[12])
col15=fits.Column(name='lnu14',format='FLOAT',array=lnu*scale[13])

cols_SFG_z2=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15])

plt.plot(lam_default,lnu*scale[0])
plt.plot(lam_default,lnu*scale[2])
plt.plot(lam_default,lnu*scale[4])
plt.plot(lam_default,lnu*scale[6])
plt.plot(lam_default,lnu*scale[8])
plt.plot(lam_default,lnu*scale[10])

plt.xscale('log')
plt.yscale('log')
plt.show()

#============================================================================
#AGN
#----------------------------------------------------------------------------

with open (seddir+'/Comprehensive_AGN1.txt','r') as f: 
    flines=f.readlines()[4:] #4 lines of header, and 4 lines as the lower-lambda values not well defined for the higher-L templates
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

lam_old=np.array(lam)
lnu_old=np.array(lnu)
if(lam_old[0] > lam_default[0]):
    lam_old=np.concatenate(([0.5],lam_old,[1031]),axis=0)
    lnu_old=np.concatenate(([0],lnu_old,[0]),axis=0)
else:
    lam_old=np.concatenate((lam_old,[1031]),axis=0)
    lnu_old=np.concatenate((lnu_old,[0]),axis=0)

f = interpolate.interp1d(lam_old,lnu_old)
lnu=f(lam_default)

col2=fits.Column(name='lnu1',format='FLOAT',array=lnu*scale[0])
col3=fits.Column(name='lnu2',format='FLOAT',array=lnu*scale[1])
col4=fits.Column(name='lnu3',format='FLOAT',array=lnu*scale[2])
col5=fits.Column(name='lnu4',format='FLOAT',array=lnu*scale[3])
col6=fits.Column(name='lnu5',format='FLOAT',array=lnu*scale[4])
col7=fits.Column(name='lnu6',format='FLOAT',array=lnu*scale[5])
col8=fits.Column(name='lnu7',format='FLOAT',array=lnu*scale[6])
col9=fits.Column(name='lnu8',format='FLOAT',array=lnu*scale[7])
col10=fits.Column(name='lnu9',format='FLOAT',array=lnu*scale[8])
col11=fits.Column(name='lnu10',format='FLOAT',array=lnu*scale[9])

#now the higher luminosity AGN template (defined from 12.1-12.8)
with open (seddir+'/Comprehensive_AGN2.txt','r') as f: 
    flines=f.readlines()[4:]
    fcount=-1
    lam,lnu=[],[]
    for fline in flines:
        tmp1,tmp2,tmp3=fline.split()
        lam.append(float(tmp1)),lnu.append(float(tmp2))

lam_old=np.array(lam)
lnu_old=np.array(lnu)
if(lam_old[0] > lam_default[0]):
    lam_old=np.concatenate(([0.5],lam_old,[1031]),axis=0)
    lnu_old=np.concatenate(([0],lnu_old,[0]),axis=0)
else:
    lam_old=np.concatenate((lam_old,[1031]),axis=0)
    lnu_old=np.concatenate((lnu_old,[0]),axis=0)

f = interpolate.interp1d(lam_old,lnu_old)
lnu=f(lam_default)

#these should be the same as above when using the same lambda
lir=0
for i in range(0,len(lam_default)-1):
    if(lam_default[i] >= 5.0):
        lir+=(lnu[i]*dnu[i])/3.86e26

print 'AGN2 template L_TIR', np.log10(lir)
scale=pow(10,lir_list-np.log10(lir))


col12=fits.Column(name='lnu11',format='FLOAT',array=lnu*scale[10])
col13=fits.Column(name='lnu12',format='FLOAT',array=lnu*scale[11])
col14=fits.Column(name='lnu13',format='FLOAT',array=lnu*scale[12])
col15=fits.Column(name='lnu14',format='FLOAT',array=lnu*scale[13])

cols_AGN_z1=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15])


#=========================================================================
# the z-1.8 set
#-------------------------------------------------------------------------

with open (seddir+'/Comprehensive_AGN3.txt','r') as f: 
    flines=f.readlines()[4:] #4 lines of header, and 4 lines as the lower-lambda values not well defined for the higher-L templates
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

lam_old=np.array(lam)
lnu_old=np.array(lnu)
if(lam_old[0] > lam_default[0]):
    lam_old=np.concatenate(([0.5],lam_old,[1031]),axis=0)
    lnu_old=np.concatenate(([0],lnu_old,[0]),axis=0)
else:
    lam_old=np.concatenate((lam_old,[1031]),axis=0)
    lnu_old=np.concatenate((lnu_old,[0]),axis=0)

f = interpolate.interp1d(lam_old,lnu_old)
lnu=f(lam_default)

lir_test=0
for i in range(0,len(lam_default)-1):
    if(lam_default[i] >= 5.0):
        lir_test+=(lnu[i]*dnu_default[i])/3.86e26

print 'AGN3 template L_TIR after interpolation', np.log10(lir_test)
scale=pow(10,lir_list-np.log10(lir_test))

col2=fits.Column(name='lnu1',format='FLOAT',array=lnu*scale[0])
col3=fits.Column(name='lnu2',format='FLOAT',array=lnu*scale[1])
col4=fits.Column(name='lnu3',format='FLOAT',array=lnu*scale[2])
col5=fits.Column(name='lnu4',format='FLOAT',array=lnu*scale[3])
col6=fits.Column(name='lnu5',format='FLOAT',array=lnu*scale[4])
col7=fits.Column(name='lnu6',format='FLOAT',array=lnu*scale[5])
col8=fits.Column(name='lnu7',format='FLOAT',array=lnu*scale[6])
col9=fits.Column(name='lnu8',format='FLOAT',array=lnu*scale[7])
col10=fits.Column(name='lnu9',format='FLOAT',array=lnu*scale[8])
col11=fits.Column(name='lnu10',format='FLOAT',array=lnu*scale[9])
col12=fits.Column(name='lnu11',format='FLOAT',array=lnu*scale[10])
col13=fits.Column(name='lnu12',format='FLOAT',array=lnu*scale[11])


#now the higher luminosity AGN template (defined from 12.1-12.8)
with open (seddir+'/Comprehensive_AGN4.txt','r') as f: 
    flines=f.readlines()[4:]
    fcount=-1
    lam,lnu=[],[]
    for fline in flines:
        tmp1,tmp2,tmp3=fline.split()
        lam.append(float(tmp1)),lnu.append(float(tmp2))


lam_old=np.array(lam)
lnu_old=np.array(lnu)
if(lam_old[0] > lam_default[0]):
    lam_old=np.concatenate(([0.5],lam_old,[1031]),axis=0)
    lnu_old=np.concatenate(([0],lnu_old,[0]),axis=0)
else:
    lam_old=np.concatenate((lam_old,[1031]),axis=0)
    lnu_old=np.concatenate((lnu_old,[0]),axis=0)

f = interpolate.interp1d(lam_old,lnu_old)
lnu=f(lam_default)

lir=0
for i in range(0,len(lam_default)-1):
    if(lam_default[i] >= 5.0):
        lir+=(lnu[i]*dnu[i])/3.86e26

print 'AGN4 template L_TIR', np.log10(lir)
scale=pow(10,lir_list-np.log10(lir))


col14=fits.Column(name='lnu13',format='FLOAT',array=lnu*scale[12])
col15=fits.Column(name='lnu14',format='FLOAT',array=lnu*scale[13])

cols_AGN_z2=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15])

#==========================================================================
# structure the different tables
#---------------------------------------------------------------------------
tbhdu10=fits.new_table(cols_SFG_z0)
tbhdu11=fits.new_table(cols_SFG_z1)
tbhdu12=fits.new_table(cols_SFG_z2)

tbhdu21=fits.new_table(cols_AGN_z1)
tbhdu22=fits.new_table(cols_AGN_z2)

tbhdu10.header.update('EXTNAME','SFG_z0',
                        'Rieke2009+SWIRE Sc')
tbhdu11.header.update('EXTNAME','SFG_z1',
                        'Kirkpatrick2015')
tbhdu12.header.update('EXTNAME','SFG_z2',
                        'Kirkpatrick2015')

tbhdu21.header.update('EXTNAME','AGN_z1',
                        'Kirkpatrick2015')
tbhdu22.header.update('EXTNAME','AGN_z2',
                        'Kirkpatrick2015')

tbhdr=tbhdu10.header
tbhdr=tbhdu11.header
tbhdr=tbhdu12.header
tbhdr=tbhdu21.header
tbhdr=tbhdu22.header

hdulist[0].header.update('Nzbins',3,'number of redshift bins')
hdulist[0].header.update('Ntypes',3,'number of SED types included')
hdulist[0].header.update('SEDTYPE1','SFG','name of 1st SED type')
hdulist[0].header.update('SEDTYPE2','AGN','name of 2nd SED type')
hdulist[0].header.update('SEDTYPE3','Comp','name of 3rd SED type')

hdulist[0].header.update('zmin0',0,'minimum of 1st redshift bin')
hdulist[0].header.update('zmax0',0.499,'maximum of 1st redshift bin')
hdulist[0].header.update('zmin1',0.5,'minimum of 2nd redshift bin')
hdulist[0].header.update('zmax1',1.2,'maximum of 2nd redshift bin')
hdulist[0].header.update('zmin2',1.201,'minimum of 3rd redshift bin')
hdulist[0].header.update('zmax2',2.5,'in code adopt the maximum z')



hdr=hdulist[0].header

hdulist.close
thdulist=fits.HDUList([hdulist[0],tbhdu10,tbhdu11,tbhdu12,tbhdu21,tbhdu22])
if(os.path.isfile(sedfile)):
    print 'Replacing existing template file...'
    os.remove(sedfile)
    thdulist.writeto(sedfile)
else:
    print 'Writing a new template file...'
    thdulist.writeto(sedfile)
