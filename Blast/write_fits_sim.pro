pro write_fits_sim,models,wave,n_params,param_sep,param_min,param_max,model_form=model_form

  if n_params() lt 6 then begin
     print,"Procedure Call Incorrect"
     print,"Syntax:"
     print,"write_fits_sim,$"
     print,"models,$ //array of all fluxes for all models"
     print,"wave,$ //array of wavelengths used to generate models (should be the size of one given model)"
     print,"n_params,$ //number of model parameters"
     print,"param_sep,$ //array of parameter step sizes"
     print,"param_min,$ //array of minima of parameters"
     print,"param_max,$ //array of maxima of parameters"
     print,"model_form=model_form,$ //string with the functional form of the model (optional)"
  endif

  models = reform(models,n_elements(models))

  sxaddpar, hdr, 'DATE', systime(),'Date of creation'
  if (n_elements(model_form) gt 0) then begin
     sxaddpar, hdr, 'MODEL',model_form,'Template Model Used'
  endif
  
  sxaddpar, hdr, 'P_NUM',n_params,'Template Model Used'

  n_els = n_elements(wave)
  wave_min = wave[0]
  wave_max = wave[n_els-1]
  wave_sep = (wave_max-wave_min)/(n_els-1)

  sxaddpar, hdr, 'WAVE_MIN',wave_min,'Domain Lower Bound'
  sxaddpar, hdr, 'WAVE_MAX',wave_max,'Domain Upper Bound'
  sxaddpar, hdr, 'WAVE_SEP',wave_sep,'Domain Step Size'
  
  for i=0,n_params-1 do begin
     p_num = strcompress(string(i),/remove_all)
     sxaddpar, hdr, 'P'+p_num+'_MAX',param_max[i],'Parameter '+p_num+' Upper Limit'
     sxaddpar, hdr, 'P'+p_num+'_MIN',param_min[i],'Parameter '+p_num+' Lower Limit'
     sxaddpar, hdr, 'P'+p_num+'_SEP',param_sep[i],'Parameter '+p_num+' Step Size'
  endfor
  
  mwrfits,models,'model.fits',hdr,/create

end
