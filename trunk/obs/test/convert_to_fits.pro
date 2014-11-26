;Nov 24 2014
; convert table with flux values to fitsfile
;Output:test.fits

pro cfits
path='/Users/andreasilva/Documents/Research_Tufts/lum_function/data/'
file='test.tbl'

;mvega and emvega are vega magnitudes with their errors for the 4 wise bands
readcol, path+file, id, ra,dec, mvega1, emvega1, mvega2, emvega2, mvega3, emvega3, mvega4, emvega4, format='A,D,D,D,D,D,D,D,D,D,D', skipline=17, /silent
n=n_elements(id)
;Convert from vega magnitudes to fluxes in mJy for the 4 WISE bands
flux1=1.0E3*10.0^((8.9-(mvega1+2.683))/2.5)
flux1u=1.0E3*10.0^((8.9-(mvega1+emvega1+2.683))/2.5)
flux1b=1.0E3*10.0^((8.9-(mvega1-emvega1+2.683))/2.5)

flux2=1.0E3*10.0^((8.9-(mvega2+2.683))/2.5)
flux2u=1.0E3*10.0^((8.9-(mvega2+emvega2+2.683))/2.5)
flux2b=1.0E3*10.0^((8.9-(mvega2-emvega2+2.683))/2.5)

flux3=1.0E3*10.0^((8.9-(mvega3+2.683))/2.5)
flux3u=1.0E3*10.0^((8.9-(mvega3+emvega3+2.683))/2.5)
flux3b=1.0E3*10.0^((8.9-(mvega3-emvega3+2.683))/2.5)

flux4=1.0E3*10.0^((8.9-(mvega4+2.683))/2.5)
flux4u=1.0E3*10.0^((8.9-(mvega4+emvega4+2.683))/2.5)
flux4b=1.0E3*10.0^((8.9-(mvega4-emvega4+2.683))/2.5)

;write into fits files
nbands=4 ;number of wise bands
wave=dblarr(nbands)
wave[0]=3.4 ;in microns
wave[1]=4.6
wave[2]=12.0
wave[3]=22.0

values=dblarr(nbands, n)
for i=0, n-1 do begin
   values(0,i)=flux1[i]
   values(1,i)=flux2[i]
   values(2,i)=flux3[i]
   values(3,i)=flux4[i]
endfor

sxaddpar, hdr, 'DATE', systime(), 'Date of creation'
sxaddpar, hdr, 'TELESCOPE WISE '
sxaddpar, hdr, 'FLUX_NUM', n, 'Number of observations included in file'
sxaddpar, hdr, 'WAVE_1',wave[0],'Wavelength corresponding to first flux'
sxaddpar, hdr, 'WAVE_2',wave[1],'Wavelength corresponding to second flux'
sxaddpar, hdr, 'WAVE_3',wave[2],'Wavelength corresponding to third flux'
sxaddpar, hdr, 'WAVE_4',wave[3],'Wavelength corresponding to fourth flux'
  

mwrfits, values, path+'test.fits', hdr, /create


end
