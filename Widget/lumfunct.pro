pro lumfunct

  phi0 = -2.2
  L0 = 10.15
  alpha = 0.47
  beta = 2.88
  p = 7
  q = 3.5

  omega_m = 0.28
  omega_l = 0.72

  set_plot,'x'
  plot,[8,13],[1e0,1e15],/ylog,/nodata
  lums = indgen(21)/4.0+8.0
  print,lums

  dz = 0.5

  for z=0.01,0.01,dz do begin
     t1 = (10^phi0)*((1+z)^p)
     t2 = (10^L0)*((1+z)^q)
     r = 10^lums/t2
     nsrcs=t1/(r^alpha+r^beta)

     ez = sqrt(omega_m*(1+z)^3+omega_l)
     dvdz = 4062*lumdist(z)^2*(!pi/180)^2/((1+z)^2*ez)
     vol = dvdz*dz
     nsrcs*=vol

     print,vol
     oplot,lums,nsrcs
  endfor

  oplot,[8,13],[1,1],linestyle=1

end
