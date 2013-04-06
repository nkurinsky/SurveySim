pro colors

  res = mrdfits(sf_templates.fits)
  help,res,/struct

  bands = [250,350,500]
  wave = res[*,0]
  inds = value_locate(wave,bands)

  set_plot'x'
  plot,[-1,1],[-1,1],/nodata,xstyle=1,ystyle=1
  
  for j=1,13 do begin
     s = res[*,j]
     for i=0.0,10.0 do begin
        b_rest = bands/(1.0+i)
        inds = value_locate(wave,b_rest)
        fluxes = s[inds]
        c1 = alog10(fluxes[2]/fluxes[0])
        c2 = alog10(fluxes[1]/fluxes[0])
        oplot,c1,c2
     end
  end
  
end
