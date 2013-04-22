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
        print,c1,c2
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

end
