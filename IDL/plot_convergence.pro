pro plot_convergence,filename

  res= mrdfits(filename,5,/silent)
  gpts = where(res.(0) gt 0)
  res = res[gpts]
  rmax = 0
  for i=0,n_elements(tag_names(res))-1 do begin
     tmax = max(res.(i))
     if (tmax gt rmax) then rmax = tmax
  endfor
  
  plot,[0,n_elements(res)],[1.0,rmax],xstyle=1,ystyle=1,title="Convergence",xtitle="Test Number",ytitle="R",/nodata
  ;oplot,[0,n_elements(gpts)],[msettings.conv_rmax,msettings.conv_rmax]
  for i=0,n_elements(tag_names(res))-1 do begin
     oplot,res.(i),linestyle=i
  endfor

end
