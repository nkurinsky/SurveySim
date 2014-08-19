pro counts
  
  res = mrdfits('output.fits',6,head,/silent)
  alpha = 0.159                ;one std deviation
  plusfrac = (1.0-alpha)
  minusfrac = alpha

  chis = res.chisq
  gpts = where(chis gt median(chis))
  res = res[gpts]

  c1=res.dnds250
  c2=res.dnds350
  c3=res.dnds500

  c1size = n_elements(res[0].dnds250)
  c2size = n_elements(res[0].dnds350)
  c3size = n_elements(res[0].dnds500)

  c1mean = make_array(c1size,value=0.0)
  c2mean = make_array(c2size,value=0.0)
  c3mean = make_array(c3size,value=0.0)
  c1plus = make_array(c1size,value=0.0)
  c2plus = make_array(c2size,value=0.0)
  c3plus = make_array(c3size,value=0.0)
  c1minus = make_array(c1size,value=0.0)
  c2minus = make_array(c2size,value=0.0)
  c3minus = make_array(c3size,value=0.0)
  
  for i=0,c1size-1 do begin
     dnds = c1[c1size-i-1,*]
     dnds = dnds[sort(dnds)]
     pi = plusfrac*n_elements(dnds)
     mi = minusfrac*n_elements(dnds)
     c1mean[i] = mean(dnds)
     c1plus[i] = dnds[pi]
     c1minus[i] = dnds[mi]
  endfor

  for i=0,c2size-1 do begin
     dnds = c2[i,*]
     dnds = dnds[sort(dnds)]
     pi = plusfrac*n_elements(dnds)
     mi = minusfrac*n_elements(dnds)
     c2mean[i] = mean(dnds)
     c2plus[i] = dnds[pi]
     c2minus[i] = dnds[mi]
  endfor

  for i=0,c3size-1 do begin
     dnds = c3[i,*]
     dnds = dnds[sort(dnds)]
     pi = plusfrac*n_elements(dnds)
     mi = minusfrac*n_elements(dnds)
     c3mean = mean(dnds)
     c3plus = dnds[pi]
     c3minus = dnds[mi]
  endfor

  set_plot,'x'

  x = indgen(n_elements(c1mean))
  plot,x,c1mean,xrange=[0,16]
  ;;oplot,x,c1plus,linestyle=1
  ;;oplot,x,c1minus,linestyle=2

end
