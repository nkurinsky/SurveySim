pro colors

  !p.thick=5
  !x.thick=5
  !y.thick=5
  !p.charthick=4
  !p.charsize=1.5

  res = mrdfits('sf_templates.fits')
  help,res,/struct

  bands = [250,350,500]
  wave = res[*,0]
  inds = value_locate(wave,bands)

  set_plot,'ps'
  device,filename='model_colors.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated
  loadct,0,/silent
  plot,[-1,1],[-0.5,0.5],/nodata,xstyle=1,ystyle=1,xtitle=textoidl("\alpha_{250}^{500}"),ytitle=textoidl("\alpha_{250}^{350}")

  loadct,39,/silent
  for j=1,14 do begin
     s = reform(res[*,j])
     for i=0.0,5.0,0.1 do begin
        b_rest = bands/(1.0+i)
        inds = value_locate(wave,b_rest)
        fluxes = s[inds]
        c1 = alog10(fluxes[2]/fluxes[0])
        c2 = alog10(fluxes[1]/fluxes[0])
        oplot,[c1],[c2],psym=2,color=i*40+20,symsize=0.5
     end
  end
  
  xyouts,-0.32,-0.38,textoidl('Redshift:')
  xyouts,-0.1,-0.38,textoidl('0.0')
  xyouts,0.09,-0.38,textoidl('1.0')
  xyouts,0.28,-0.38,textoidl('2.0')
  xyouts,0.47,-0.38,textoidl('3.0')
  xyouts,0.66,-0.38,textoidl('4.0')
  xyouts,0.85,-0.38,textoidl('5.0')
  for i=0.0,5.0,0.1 do begin
     ind = (-0.1+i/5)
     oplot,[ind,ind+0.01,ind,ind+0.01,ind,ind+0.01],[-0.4,-0.4,-0.395,-0.395,-0.39,-0.39],psym=2,color=i*40+20,symsize=0.5
  endfor

  device,/close
  device,filename="models.eps",xsize=12,ysize=9,/inches,/times,/color,/encapsulated
  loadct,0,/silent

  plot,[wave[0],wave[n_elements(wave)-1]],[1e19,1e28],/xlog,/ylog,/nodata,xstyle=1,ystyle=1,xtitle=textoidl("Wavelength [\mu m]"),ytitle=textoidl("Luminosity [W/\lambda]")

  for j=1,14 do begin
     s = reform(res[*,j])
     oplot,wave,s,linestyle=0
  endfor

  oplot,[250,250],[1e19,1e28],linestyle=1
  oplot,[350,350],[1e19,1e28],linestyle=1
  oplot,[500,500],[1e19,1e28],linestyle=1

  device,/close

  device,filename='model_brightness.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated
  
  plot,[5e1,5e3],[1e19,1e28],/xlog,/ylog,/nodata,xstyle=1,ystyle=1,xtitle=textoidl("Wavelength [\mu m]"),ytitle=textoidl("Luminosity [W/\lambda]")

  loadct,39,/silent
  for j=1,1 do begin
     s = reform(res[*,j])
     for i=0.0,5.0 do begin
        b_obs = wave*(1.0+i)
        oplot,b_obs,s,color=i*40+20
     endfor
  endfor

  s = reform(res[*,14])
  for i=0.0,5.0 do begin
     b_obs = wave*(1.0+i)
     oplot,b_obs,s,color=i*40+20
  endfor

  for i=0,2 do begin
     oplot,[bands[i],bands[i]],[1e19,1e28],linestyle=1
  endfor	

  xyouts,6e1,10^(20.7),textoidl('Redshift:')
  xyouts,6e1,10^(20.4),textoidl('0.0')
  xyouts,2e2,10^(20.4),textoidl('5.0')
  for i=0.0,5.0,0.1 do begin
     ind = (1.8+i/10)
     oplot,[10^ind,10^(ind+0.005),10^ind,10^(ind+0.005)],[10^20.0,10^20.0,10^20.1,10^20.1],psym=2,color=i*40+20,symsize=1.0
  endfor

  device,/close

end
