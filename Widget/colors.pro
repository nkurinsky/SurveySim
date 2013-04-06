pro colors

  res = mrdfits('sf_templates.fits')
  help,res,/struct

  bands = [250,350,500]
  wave = res[*,0]
  inds = value_locate(wave,bands)

  set_plot,'x'
  window,xsize=800,ysize=600
  plot,[-1,1],[-0.5,0.5],/nodata,xstyle=1,ystyle=1

  col1 = fltarr(14*11)
  col2 = col1
  
  for j=1,14 do begin
     s = reform(res[*,j])
     for i=0.0,10.0 do begin
        b_rest = bands/(1.0+i)
        inds = value_locate(wave,b_rest)
        print,inds
        fluxes = s[inds]
        c1 = alog10(fluxes[2]/fluxes[0])
        c2 = alog10(fluxes[1]/fluxes[0])
        print,c1,c2
        oplot,[c1],[c2],psym=2
     end
  end
  
end
