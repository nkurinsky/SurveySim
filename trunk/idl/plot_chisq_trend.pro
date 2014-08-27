pro plot_chisq_trend,filename
  
  res = mrdfits(filename,4,/silent)
  alltags = tag_names(res)
  tags = alltags[where(strmatch(alltags,'*0',/FOLD_CASE) eq 1)]
  
  for i=0,n_elements(tags)-1 do begin
     tags[i] = strmid(tags[i],0,strlen(tags[i])-1)
  endfor
  tags = tags[where((tags ne 'CHISQ') and (tags ne 'ACPT'))]
  print,'Fitted Variables: ',tags

  dim = n_elements(tags)
  nbins = 20.0
  
  cpts = where(strmatch(alltags,'CHISQ*',/FOLD_CASE) eq 1)
  for ci=0,n_elements(cpts)-1 do begin
     if(ci eq 0) then begin
        chis = res.(cpts[ci])
     endif else begin
        chis = [chis,res.(cpts[ci])]
     endelse
  endfor
  
  apts = where(strmatch(alltags,'ACPT*',/FOLD_CASE) eq 1)
  accept = make_array(n_elements(apts),n_elements(res.(cpts[0])),value=0.0)
  for ai=0,n_elements(apts)-1 do begin
     if(ai eq 0) then begin
        accept = res.(apts[ai])
     endif else begin
        accept = [accept,res.(apts[ai])]
     endelse
  endfor

  gpts = where(accept eq 1.0)
  crange=[min(chis[gpts]),max(chis[gpts])]
  n = n_elements(res)
  
  loadct,1,/silent
  plot,[0,n],crange,/ylog,xstyle=1,/nodata,xtitle='Run Number',ytitle=textoidl("\chi^2"),title="Temporal Likelihood Trends"
  loadct,39,/silent
  
  apts = where(strmatch(alltags,'ACPT*',/FOLD_CASE) eq 1)
  cpts = where(strmatch(alltags,'CHISQ*',/FOLD_CASE) eq 1)
  chainnum = n_elements(apts)
  dcolor = 200/(chainnum-1)
  xchis = indgen(n)
  for i=0,chainnum-1 do begin
     gpts = where(res.(apts[i]) eq 1.0)
     yplot = res.(cpts[i])
     oplot,xchis[gpts],yplot[gpts],color=40+i*dcolor,psym=3
  endfor

end
