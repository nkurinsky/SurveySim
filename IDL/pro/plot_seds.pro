pro plot_seds,savefile="params.save"

  COMMON simulation_com
  
  set_plot,'x'
  device,decomposed=0
  plot_settings,plot_type='x'

  if(n_elements(parameters) eq 0) then begin
     print,"Using plot_seds in standalone mode, looking for "+savefile
     parameters=load_parameters(savefile)
  endif else begin
     wset,info.wind_id3
  endelse
  
  sedfile = parameters.files.sedfile

  if file_test(sedfile) then begin
     templ=mrdfits(sedfile,/silent) 
     loadct,1,/silent
     plot,templ[*,0],templ[*,1],/xlog,/ylog,yrange=[1.d20,1.d28],ystyle=1,xtitle=TeXtoIDL('\lambda [\mum]'),ytitle=TeXtoIDL('L_{\nu} [W/Hz]')
     for ipl=1,13 do oplot,templ[*,0],templ[*,ipl+1]
  endif else begin
     print,'File '+sedfile+' does not exist, skipping plot'
  endelse
  
end
