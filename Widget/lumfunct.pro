pro lumfunct

  !p.thick=5
  !x.thick=5
  !y.thick=5
  !p.charthick=4
  !p.charsize=1.5

  phi0 = -1.8
  L0 = 10.15
  alpha = 0.47
  beta = 2.88
  p = -6.7
  q = 3.5

  omega_m = 0.28
  omega_l = 0.72

  set_plot,'ps'
  device,filename='lumfunct.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated

  plot,[8,13],[1e-5,1e5],/ylog,/nodata,xtitle=textoidl('log_{10}(L_{fir}) [L_{sun}]'),ytitle=textoidl('log_{10}(N(L_{fir})/\Omega)'),ystyle=1
  lums = indgen(21)/4.0+8.0
  print,lums

  dz = 0.25

  loadct,39,/silent
  for z=dz,5.0,dz do begin
     if z le 2.0 then begin
        t1 = (10^phi0)*((1+z)^p)
        t2 = (10^L0)*((1+z)^q)
     endif else begin
        t1 = (10^phi0)*((3.0)^p)
        t2 = (10^L0)*((3.0)^q)
     endelse

     r = 10^lums/t2
     nsrcs=t1/(r^alpha+r^beta)

     ez = sqrt(omega_m*(1+z)^3+omega_l)
     dvdz = 4062*lumdist(z)^2*(!pi/180)^2/((1+z)^2*ez)
     vol = dvdz*dz
     nsrcs*=vol

     print,vol
     oplot,lums,nsrcs,color=z*40+20
  endfor

  for i=0.0,5.0,0.1 do begin
     ind = (9+2*i/5)
     oplot,[ind,ind+0.02,ind,ind+0.02,ind,ind+0.02],[10^(-3),10^(-3),10^(-2.9),10^(-2.9),10^(-2.8),10^(-2.8)],psym=2,color=i*40+20,symsize=0.75
  endfor

  loadct,0,/silent
  oplot,[8,13],[1,1],linestyle=1

  xyouts,9.0,10^(-2.2),textoidl('Redshift')
  xyouts,9.0,10^(-2.6),textoidl('0')
  xyouts,10.9,10^(-2.6),textoidl('5')

  device,/close

end
