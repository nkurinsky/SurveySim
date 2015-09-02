#!/usr/bin/env python

import os
import numpy as np
from astropy.io import ascii
from astropy.io import fits

#name of the output fits file (include full path)
outfits='/Users/annie/students/noah_kurinsky/SurveySim/trunk/obs/wise.fits'

#ascii file contains the magnitudes (and their errors) for the 4 wise bands. 
#given in  Vega magnitudes
file='/Users/annie/students/noah_kurinsky/SurveySim/trunk/obs/test/test.tbl'
data=ascii.read(file);

#Convert from vega magnitudes to fluxes in mJy for the 4 WISE bands
flux1=1.0E3*pow(10.0,(8.9-(np.array(data['w1mpro'])+2.699))/2.5);
flux2=1.0E3*pow(10.0,(8.9-(np.array(data['w2mpro'])+3.339))/2.5);
flux3=1.0E3*pow(10.0,(8.9-(np.array(data['w3mpro'])+5.174))/2.5);
flux4=1.0E3*pow(10.0,(8.9-(np.array(data['w4mpro'])+6.620))/2.5);

eflux1=flux1*abs((1-pow(10,(-0.4*np.array(data['w1sigmpro'])))))
eflux2=flux2*abs((1-pow(10,(-0.4*np.array(data['w2sigmpro'])))))
eflux3=flux3*abs((1-pow(10,(-0.4*np.array(data['w3sigmpro'])))))
eflux4=flux4*abs((1-pow(10,(-0.4*np.array(data['w4sigmpro'])))))

col1 = fits.Column(name='F3.4', format='D', array=flux1,unit='mJy');
col2 = fits.Column(name='EF3.4', format='D', array=eflux1,unit='mJy');
col3 = fits.Column(name='F4.6', format='D', array=flux2,unit='mJy');
col4 = fits.Column(name='EF4.6', format='D', array=eflux2,unit='mJy');
col5 = fits.Column(name='F12', format='D', array=flux3,unit='mJy');
col6 = fits.Column(name='EF12', format='D', array=eflux3,unit='mJy');
col7 = fits.Column(name='F22', format='D', array=flux4,unit='mJy');
col8 = fits.Column(name='EF22', format='D', array=eflux4,unit='mJy');

cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8]);
tbhdu = fits.new_table(cols);

#create a new fits file with a binary table in its 1st extension
hdu = fits.PrimaryHDU()
#make any necessary updates to the primary HDU header
prihdr=hdu.header
prihdr.set('FHDU',1,'location of data table')
#make any necessary updates to the table extension header
thdr=tbhdu.header
thdr.set('NBANDS',4,'number of bands')
thdr.set('BAND1',3.4,'in microns')
thdr.set('BAND2',4.6, 'in microns')
thdr.set('BAND3',12.0,'in microns')
thdr.set('BAND4',22.0,'in microns')

thdr.set('F1MIN',min(flux1),'Minimum flux1, these data')
thdr.set('F2MIN',min(flux2),'Minimum flux2, these data')
thdr.set('F3MIN',min(flux3),'Minimum flux3, these data')
thdr.set('F4MIN',min(flux4),'Minimum flux4, these data')
thdr.set('F1ERR',np.mean(eflux1),'rms1, these data')
thdr.set('F2ERR',np.mean(eflux2),'rms2, these data ')
thdr.set('F3ERR',np.mean(eflux3),'rms3, these data')
thdr.set('F4ERR',np.median(eflux4),'rms4, these data')

#don't really want these to be hardwired here. should remove the bits of the cpp code that require them (obs_lib.cpp)
thdr.set('F1COL','F3.4','Name of the flux1 column')
thdr.set('F2COL','F4.6','Name of the flux2 column')
thdr.set('F3COL','F12','Name of the flux3 column')
thdr.set('F4COL','F22','Name of the flux4 column')
thdr.set('EF1COL','EF3.4','Name of the eflux1 column')
thdr.set('EF2COL','EF4.6','Name of the eflux2 column')
thdr.set('EF3COL','EF12','Name of the eflux3 column')
thdr.set('EF4COL','EF22','Name of the eflux4 column')
thdr.set('F1FILT','W1_3.4')
thdr.set('F2FILT','W2_4.6')
thdr.set('F3FILT','W3_12')
thdr.set('F4FILT','W4_22') 

thdulist = fits.HDUList([hdu, tbhdu])
#check if output file already exists and replace with new one if so
if(os.path.isfile(outfits)):
   print 'Replacing existing file....'
   os.remove(outfits)
   thdulist.writeto(outfits)
else:
   thdulist.writeto(outfits)






