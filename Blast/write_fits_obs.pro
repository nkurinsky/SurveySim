pro write_fits_obs,f1,f2,f3,wave,observatory=observatory,filename=filename

  if (n_elements(wave) ne 3) then begin
     print,'"wave" must be a vector of size=3'
     return
  endif
  
  if n_params() lt 4 then begin
     print,"Procedure Call Incorrect"
     print,"Syntax:"
     print,"write_fits_obs,$"
     print,"f1,f2,f3,$ //array of fluxes for each of the three bands"
     print,"wave,$ //array of wavelengths used to generate models (should be the size of one given model)"
     print,"observatory=observatory$ //string with the name of observatory"
     return
  endif

  values = dblarr(n_elements(wave),n_elements(f1))
  
  for i=0,n_elements(f1)-1 do begin
     values(0,i) = f1[i]
     values(1,i) = f2[i]
     values(2,i) = f3[i]
  endfor

  sxaddpar, hdr, 'DATE', systime(),'Date of creation'
  if (n_elements(observatory) gt 0) then begin
     sxaddpar, hdr, 'TELESCOP',observatory,'Observatory'
  endif
  
  sxaddpar, hdr, 'FLUX_NUM',n_elements(f1),'Number of observations included in file'
  sxaddpar, hdr, 'WAVE_1',wave[0],'Wavelength corresponding to first flux'
  sxaddpar, hdr, 'WAVE_2',wave[1],'Wavelength corresponding to second flux'
  sxaddpar, hdr, 'WAVE_3',wave[2],'Wavelength corresponding to third flux'
  
  if(n_elements(filename) eq 0) then begin
     filename = 'observation.fits'
  endif

  for i=0,5 do begin
     print,values[i]
  endfor

  mwrfits,values,filename,hdr,/create

end
