;Makes SEDs of some of the extreme sources from observation
;Current as of 6/22

pro sed_tests

  !p.multi = [0,2,2]
  restore,'save_files/matched_fluxes.save'
  c1 = acolor(f250,f350,250.d0,350.d0)
  c2 = acolor(f350,f500,350.d0,500.d0)
  bands = [250.d0,350.d0,500.d0]
  true_bands = bands*1d-6
  types = ['dense1','dense2','scarce1','scarce2']
  
  d1 = where((c1 gt -0.5) and (c1 lt 0) and (c2 lt -1) and (c2 gt -1.5))
  d2 = where((c1 gt -1) and (c1 lt -0.5) and (c2 lt -0.5) and (c2 gt -1))
  s1 = where((c1 gt 1) and (c2 lt -2.5))
  s2 = where((c1 lt -2) and (c2 gt 1))
  inds = [d1[0],d2[0],s1[0],s2[0]]

  set_plot,'ps'
  device,filename='plots/seds.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  for i=0,3 do begin
     fluxes = [f250[inds[i]],f350[inds[i]],f500[inds[i]]]
     plot,bands,fluxes,title=types[i],xrange=[100,1000],/xlog,xstyle=1,/ylog,ystyle=1,yrange=[1e-10,1e15]
     params = [100.d0,2.0d]
     res = curvefit(true_bands,fluxes,replicate(1.0,n_elements(fluxes)),params,errs,function_name='bbfunct')
  endfor
  device,/close
  
end
