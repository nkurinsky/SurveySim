#!/usr/bin/env python

#Python prelims
import os
import sys

#codedir=os.getcwd()+'/trunk/';
pydir=os.getcwd();
#codedir = pydir.replace(' ', '')[:-6] #remove last 6 characters (i.e. Python)
seddir=os.getcwd()+'/trunk/templates/TSV'

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

sedfile=os.getcwd()+'/trunk/templates/default_templates.fits'

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

#want to extend to 0.2um 
lam_short=[0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.10,4.20,4.30,4.40,4.50,4.60,4.70,5.0,5.1,5.2]
lam_short=np.array(lam_short)
lam_rieke=np.array(lam_rieke)
lam_default=np.concatenate((lam_short,lam_rieke),axis=0)

dnu_default=(np.abs(lam_default-np.roll(lam_default,1))/lam_default)/lam_default 
dnu_default=(c*1e6)*dnu_default

#read in Sc spiral galaxy template from SWIRE template library in order to get the stellar bump and switch to the above at 0.5-4um
#use the Sc template for the low-L cases (logL<11.7)
with open ('/Users/annie/code/idl/auxdata/swire/Sc_template_norm.sed','r') as f:
    flines=f.readlines()
    fcount=-1
    lam_swire,flam_swire=[],[]
    for fline in flines:
        tmp0,tmp1=fline.split()
        lam_swire.append(float(tmp0)),flam_swire.append(float(tmp1))

#use the IRAS20551 starburst template for the higher-L cases
with open ('/Users/annie/code/idl/auxdata/swire/I20551_template_norm.sed','r') as f:
    flines=f.readlines()
    fcount=-1
    lam_swire_ulirg,flam_swire_ulirg=[],[]
    for fline in flines:
        tmp0,tmp1=fline.split()
        lam_swire_ulirg.append(float(tmp0)),flam_swire_ulirg.append(float(tmp1))

lam_swire=np.array(lam_swire)/1e4 #convert to microns
lnu_swire=np.array(flam_swire)*lam_swire**2. #convert to fnu but overall normalization not set yet!

lam_swire_ulirg=np.array(lam_swire_ulirg)/1e4 #convert to microns
lnu_swire_ulirg=np.array(flam_swire_ulirg)*lam_swire_ulirg**2. #convert to fnu but overall normalization not set yet!

f_tmp=interpolate.interp1d(lam_swire_ulirg,lnu_swire_ulirg)
lnu_swire_ulirg_interp=f_tmp(lam_swire)

f = interpolate.interp1d(lam_swire,lnu_swire)
f_ulirg=interpolate.interp1d(lam_swire,lnu_swire_ulirg_interp)
#print lam_swire[0],lam_swire[1196]
lnu_tmp=f(lam_short)
lnu_tmp_ulirg=f_ulirg(lam_short)
istitch=45
#print lam_short[istitch],lnu_tmp[istitch],lam_rieke[0],lnu1[0]
tmp=lnu_tmp*(lnu1[0]/lnu_tmp[istitch])
lnu1_default=np.concatenate((tmp,lnu1),axis=0)
tmp=lnu_tmp*(lnu2[0]/lnu_tmp[istitch])
lnu2_default=np.concatenate((tmp,lnu2),axis=0)
tmp=lnu_tmp*(lnu3[0]/lnu_tmp[istitch])
lnu3_default=np.concatenate((tmp,lnu3),axis=0)
tmp=lnu_tmp*(lnu4[0]/lnu_tmp[istitch])
lnu4_default=np.concatenate((tmp,lnu4),axis=0)
tmp=lnu_tmp*(lnu5[0]/lnu_tmp[istitch])
lnu5_default=np.concatenate((tmp,lnu5),axis=0)
tmp=lnu_tmp*(lnu6[0]/lnu_tmp[istitch])
lnu6_default=np.concatenate((tmp,lnu6),axis=0)
tmp=lnu_tmp*(lnu7[0]/lnu_tmp[istitch])
lnu7_default=np.concatenate((tmp,lnu7),axis=0)

#from here on use the ulirg-based near-IR SED
tmp=lnu_tmp_ulirg*(lnu8[0]/lnu_tmp_ulirg[istitch])
lnu8_default=np.concatenate((tmp,lnu8),axis=0)
tmp=lnu_tmp_ulirg*(lnu9[0]/lnu_tmp_ulirg[istitch])
lnu9_default=np.concatenate((tmp,lnu9),axis=0)
tmp=lnu_tmp_ulirg*(lnu10[0]/lnu_tmp_ulirg[istitch])
lnu10_default=np.concatenate((tmp,lnu10),axis=0)
tmp=lnu_tmp_ulirg*(lnu11[0]/lnu_tmp_ulirg[istitch])
lnu11_default=np.concatenate((tmp,lnu11),axis=0)
tmp=lnu_tmp_ulirg*(lnu12[0]/lnu_tmp_ulirg[istitch])
lnu12_default=np.concatenate((tmp,lnu12),axis=0)
tmp=lnu_tmp_ulirg*(lnu13[0]/lnu_tmp_ulirg[istitch])
lnu13_default=np.concatenate((tmp,lnu13),axis=0)
tmp=lnu_tmp_ulirg*(lnu14[0]/lnu_tmp_ulirg[istitch])
lnu14_default=np.concatenate((tmp,lnu14),axis=0)

dnu=(np.abs(lam_default-np.roll(lam_default,1))/lam_default)/lam_default 
dnu=(c*1e6)*dnu
lir1=0
lir2=0
lir3=0
lir4=0
lir5=0
lir6=0
lir7=0
lir8=0
lir9=0
lir10=0
lir11=0
lir12=0
lir13=0
lir14=0

for i in range(0,len(lam_default)-1):
    if(lam_default[i] >= 5.0):
        lir1+=(lnu1_default[i]*dnu[i])/3.86e26
        lir2+=(lnu2_default[i]*dnu[i])/3.86e26
        lir3+=(lnu3_default[i]*dnu[i])/3.86e26
        lir4+=(lnu4_default[i]*dnu[i])/3.86e26
        lir5+=(lnu5_default[i]*dnu[i])/3.86e26
        lir6+=(lnu6_default[i]*dnu[i])/3.86e26
        lir7+=(lnu7_default[i]*dnu[i])/3.86e26
        lir8+=(lnu8_default[i]*dnu[i])/3.86e26
        lir9+=(lnu9_default[i]*dnu[i])/3.86e26
        lir10+=(lnu10_default[i]*dnu[i])/3.86e26
        lir11+=(lnu11_default[i]*dnu[i])/3.86e26
        lir12+=(lnu12_default[i]*dnu[i])/3.86e26
        lir13+=(lnu13_default[i]*dnu[i])/3.86e26
        lir14+=(lnu14_default[i]*dnu[i])/3.86e26


print 'SFG0 template L_TIR', np.log10(lir1)
scale1=pow(10,lir_list[0]-np.log10(lir1))
scale2=pow(10,lir_list[1]-np.log10(lir2))
scale3=pow(10,lir_list[2]-np.log10(lir3))
scale4=pow(10,lir_list[3]-np.log10(lir4))
scale5=pow(10,lir_list[4]-np.log10(lir5))
scale6=pow(10,lir_list[5]-np.log10(lir6))
scale7=pow(10,lir_list[6]-np.log10(lir7))
scale8=pow(10,lir_list[7]-np.log10(lir8))
scale9=pow(10,lir_list[8]-np.log10(lir9))
scale10=pow(10,lir_list[9]-np.log10(lir10))
scale11=pow(10,lir_list[10]-np.log10(lir11))
scale12=pow(10,lir_list[11]-np.log10(lir12))
scale13=pow(10,lir_list[12]-np.log10(lir13))
scale14=pow(10,lir_list[13]-np.log10(lir14))

col1=fits.Column(name='lambda',format='FLOAT',array=lam_default)
col2=fits.Column(name='lnu1',format='FLOAT',array=lnu1_default*scale1)
col3=fits.Column(name='lnu2',format='FLOAT',array=lnu2_default*scale2)
col4=fits.Column(name='lnu3',format='FLOAT',array=lnu3_default*scale3)
col5=fits.Column(name='lnu4',format='FLOAT',array=lnu4_default*scale4)
col6=fits.Column(name='lnu5',format='FLOAT',array=lnu5_default*scale5)
col7=fits.Column(name='lnu6',format='FLOAT',array=lnu6_default*scale6)
col8=fits.Column(name='lnu7',format='FLOAT',array=lnu7_default*scale7)
col9=fits.Column(name='lnu8',format='FLOAT',array=lnu8_default*scale8)
col10=fits.Column(name='lnu9',format='FLOAT',array=lnu9_default*scale9)
col11=fits.Column(name='lnu10',format='FLOAT',array=lnu10_default*scale10)
col12=fits.Column(name='lnu11',format='FLOAT',array=lnu11_default*scale11)
col13=fits.Column(name='lnu12',format='FLOAT',array=lnu12_default*scale12)
col14=fits.Column(name='lnu13',format='FLOAT',array=lnu13_default*scale13)
col15=fits.Column(name='lnu14',format='FLOAT',array=lnu14_default*scale14)

#plt.plot(lam_default,lnu4_default)
#plt.xscale('log')
#plt.yscale('log')
#plt.show()
#print len(lam_default),len(lnu1_default)

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
    lam_old=np.concatenate(([lam_default[0]],lam_old,[1031]),axis=0)
    lnu_old=np.concatenate(([0],lnu_old,[0]),axis=0)
else:
    lam_old=np.concatenate((lam_old,[1031]),axis=0)
    lnu_old=np.concatenate((lnu_old,[0]),axis=0)

f = interpolate.interp1d(lam_old,lnu_old)
#print lam_old[0],lam_default[0]
lnu=f(lam_default)

dnu=(np.abs(lam_default-np.roll(lam_default,1))/lam_default)/lam_default 
dnu=(c*1e6)*dnu
lir=0
for i in range(0,len(lam_default)-1):
    if(lam_default[i] >= 5.0):
        lir+=(lnu[i]*scale[0]*dnu[i])/3.86e26

print 'SFG1 template L_TIR (after interp)', np.log10(lir)

#this is if we scale the existing SFG1 template from Kirkpatrick+2015 down to logL=9.75
#col2=fits.Column(name='lnu1',format='FLOAT',array=lnu*scale[0])
#col3=fits.Column(name='lnu2',format='FLOAT',array=lnu*scale[1])
#col4=fits.Column(name='lnu3',format='FLOAT',array=lnu*scale[2])
#col5=fits.Column(name='lnu4',format='FLOAT',array=lnu*scale[3])
#col6=fits.Column(name='lnu5',format='FLOAT',array=lnu*scale[4])
#col7=fits.Column(name='lnu6',format='FLOAT',array=lnu*scale[5])
#col8=fits.Column(name='lnu7',format='FLOAT',array=lnu*scale[6])

#this is if we assume that the z=0 templates for low-L galaxies are constant
col2=fits.Column(name='lnu1',format='FLOAT',array=lnu1_default*scale1)
col3=fits.Column(name='lnu2',format='FLOAT',array=lnu2_default*scale2)
col4=fits.Column(name='lnu3',format='FLOAT',array=lnu3_default*scale3)
col5=fits.Column(name='lnu4',format='FLOAT',array=lnu4_default*scale4)
col6=fits.Column(name='lnu5',format='FLOAT',array=lnu5_default*scale5)
col7=fits.Column(name='lnu6',format='FLOAT',array=lnu6_default*scale6)
col8=fits.Column(name='lnu7',format='FLOAT',array=lnu7_default*scale7)

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
    lam_old=np.concatenate(([lam_default[0]],lam_old,[1031]),axis=0)
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
    lam_old=np.concatenate(([lam_default[0]],lam_old,[1031]),axis=0)
    lnu_old=np.concatenate(([0],lnu_old,[0]),axis=0)
else:
    lam_old=np.concatenate((lam_old,[1031]),axis=0)
    lnu_old=np.concatenate((lnu_old,[0]),axis=0)

f = interpolate.interp1d(lam_old,lnu_old)
lnu=f(lam_default)

#scaling from Kirkpatrick+2015 SFGz2
#col2=fits.Column(name='lnu1',format='FLOAT',array=lnu*scale[0])
#col3=fits.Column(name='lnu2',format='FLOAT',array=lnu*scale[1])
#col4=fits.Column(name='lnu3',format='FLOAT',array=lnu*scale[2])
#col5=fits.Column(name='lnu4',format='FLOAT',array=lnu*scale[3])
#col6=fits.Column(name='lnu5',format='FLOAT',array=lnu*scale[4])
#col7=fits.Column(name='lnu6',format='FLOAT',array=lnu*scale[5])
#col8=fits.Column(name='lnu7',format='FLOAT',array=lnu*scale[6])

#this is if we assume that the z=0 templates for low-L galaxies are constant
col2=fits.Column(name='lnu1',format='FLOAT',array=lnu1_default*scale1)
col3=fits.Column(name='lnu2',format='FLOAT',array=lnu2_default*scale2)
col4=fits.Column(name='lnu3',format='FLOAT',array=lnu3_default*scale3)
col5=fits.Column(name='lnu4',format='FLOAT',array=lnu4_default*scale4)
col6=fits.Column(name='lnu5',format='FLOAT',array=lnu5_default*scale5)
col7=fits.Column(name='lnu6',format='FLOAT',array=lnu6_default*scale6)
col8=fits.Column(name='lnu7',format='FLOAT',array=lnu7_default*scale7)

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

#for z=0 use Mullaney et al template. 
mullaney_filename='/Users/annie/code/idl/auxdata/mullaney_table3.tex'
with open (mullaney_filename,'r') as f:
    flines=f.readlines()[6:]
    agn_z0_lam,agn_z0_lnu=[],[]
    for fline in flines:
        tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11=fline.split()
        agn_z0_lam.append(float(tmp1)),agn_z0_lnu.append(float(tmp2))


f = interpolate.interp1d([lam_default[0],agn_z0_lam],[0,agn_z0_lnu])
lnu=f(lam_default)
#these should be the same as above when using the same lambda
#dnu=(np.abs(lam-np.roll(lam,1))/lam)/lam 
#dnu=(c*1e6)*dnu
lir=0
for i in range(0,len(lam_default)-1):
    if(lam_default[i] >= 5.0):
        lir+=(lnu[i]*dnu[i])/3.86e26

print 'AGN0 template L_TIR', np.log10(lir)
scale=pow(10,lir_list-np.log10(lir))

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


cols_AGN_z0=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15])

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
    lam_old=np.concatenate(([lam_default[0]],lam_old,[1031]),axis=0)
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
    lam_old=np.concatenate(([lam_default[0]],lam_old,[1031]),axis=0)
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
    lam_old=np.concatenate(([lam_default[0]],lam_old,[1031]),axis=0)
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
    lam_old=np.concatenate(([lam_default[0]],lam_old,[1031]),axis=0)
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
print len(lam_default),len(lnu)
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

tbhdu10.header.update('EXTNAME','AGN_z0',
                        'Mullaney2011')
tbhdu21.header.update('EXTNAME','AGN_z1',
                        'Kirkpatrick2015')
tbhdu22.header.update('EXTNAME','AGN_z2',
                        'Kirkpatrick2015')

tbhdr=tbhdu10.header
tbhdr=tbhdu11.header
tbhdr=tbhdu12.header
tbhdr=tbhdu10.header
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
hdulist[0].header.update('SCALE',-6,'wavelength scale relative to m')

hdulist[0].header.update('L1',lir_list[0],'luminosity array')
hdulist[0].header.update('L2',lir_list[1],'luminosity array')
hdulist[0].header.update('L3',lir_list[2],'luminosity array')
hdulist[0].header.update('L4',lir_list[3],'luminosity array')
hdulist[0].header.update('L5',lir_list[4],'luminosity array')
hdulist[0].header.update('L6',lir_list[5],'luminosity array')
hdulist[0].header.update('L7',lir_list[6],'luminosity array')
hdulist[0].header.update('L8',lir_list[7],'luminosity array')
hdulist[0].header.update('L9',lir_list[8],'luminosity array')
hdulist[0].header.update('L10',lir_list[9],'luminosity array')
hdulist[0].header.update('L11',lir_list[10],'luminosity array')
hdulist[0].header.update('L12',lir_list[11],'luminosity array')
hdulist[0].header.update('L13',lir_list[12],'luminosity array')
hdulist[0].header.update('L14',lir_list[13],'luminosity array')


hdr=hdulist[0].header

hdulist.close
thdulist=fits.HDUList([hdulist[0],tbhdu10,tbhdu11,tbhdu12,tbhdu10,tbhdu21,tbhdu22])
if(os.path.isfile(sedfile)):
    print 'Replacing existing template file...'
    os.remove(sedfile)
    thdulist.writeto(sedfile)
else:
    print 'Writing a new template file...'
    thdulist.writeto(sedfile)
