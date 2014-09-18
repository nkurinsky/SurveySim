function compute_counts,filename

  count_dists = mrdfits(filename,6,head,/silent)
  
  alpha = 0.159                 ;one std deviation
  plusfrac = (1.0-alpha)
  minusfrac = alpha

  chis = count_dists.chisq
  ;gpts = where(chis lt median(chis))
  ;count_dists = count_dists[gpts]
  c1=count_dists.(1)
  c2=count_dists.(2)
  c3=count_dists.(3)

  c1size = n_elements(c1[*,0])
  c2size = n_elements(c2[*,0])
  c3size = n_elements(c3[*,0])

  counts1 = {Counts1,mean:make_array(c1size,value=0.0),plus:make_array(c1size,value=0.0),minus:make_array(c1size,value=0.0)}
  counts2 = {Counts2,mean:make_array(c2size,value=0.0),plus:make_array(c2size,value=0.0),minus:make_array(c2size,value=0.0)}
  counts3 = {Counts3,mean:make_array(c3size,value=0.0),plus:make_array(c3size,value=0.0),minus:make_array(c3size,value=0.0)}

  for i=0,c1size-1 do begin
     dnds = c1[c1size-i-1,*]
     dnds = dnds[sort(dnds)]
     dnds = dnds[where(dnds gt 0)]
     pi = plusfrac*n_elements(dnds)
     mi = minusfrac*n_elements(dnds)
     counts1.mean[i] = mean(dnds)
     counts1.plus[i] = dnds[pi]
     counts1.minus[i] = dnds[mi]
  endfor

  for i=0,c2size-1 do begin
     dnds = c2[i,*]
     dnds = dnds[sort(dnds)]
     dnds = dnds[where(dnds gt 0)]
     pi = plusfrac*n_elements(dnds)
     mi = minusfrac*n_elements(dnds)
     counts2.mean[i] = mean(dnds)
     counts2.plus[i] = dnds[pi]
     counts2.minus[i] = dnds[mi]
  endfor

  for i=0,c3size-1 do begin
     dnds = c3[i,*]
     dnds = dnds[sort(dnds)]
     dnds = dnds[where(dnds gt 0)]
     pi = plusfrac*n_elements(dnds)
     mi = minusfrac*n_elements(dnds)
     counts3.mean[i] = mean(dnds)
     counts3.plus[i] = dnds[pi]
     counts3.minus[i] = dnds[mi]
  endfor

  counts = {NumberCounts,band1:counts1,band2:counts2,band3:counts3}
  
  return,counts
end
