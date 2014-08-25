pro plot_simulation_hist,filename,extension
  
  head = headfits(filename)
  res = mrdfits(filename,extension,/silent)
  
  xysize = fxpar(head,'DIM')
  hist_min = fxpar(head,'H_MIN')
  hist_max = fxpar(head,'H_MAX')
  binsize = fxpar(head,'BINSIZE')
  
  a = findgen(xysize+1)
  a *= binsize
  a += hist_min
  
  pos = where(res gt 0)
  neg = where(res lt 0)
  resmod = res
  resobs = abs(res)
  if(neg[0] ne -1) then resmod[neg] = 0
  if(pos[0] ne -1) then resobs[pos] = 0
  plot,[hist_min,hist_max],[hist_min,hist_max],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}')
  hb = 0.5

  if(extension ne 0) then begin
     color = 260*res/(max(res)*1.2)+30
     color[where(res eq 0)] = 0
     loadct,39,/silent
     for i=0,xysize-1 do begin
        for j=0,xysize-1 do begin
           xfill = [a[i],a[i],a[i+1],a[i+1]]
           yfill = [a[j],a[j+1],a[j+1],a[j]]
           if(res(i,j) gt 0) then begin
              polyfill,xfill,yfill,color=color(i,j)
           endif
        endfor
     endfor
  endif else begin
     for resi=0,1 do begin
        case resi of
           0:begin
              wres = resmod
              loadct,3,/silent
           end
           1:begin
              wres = resobs
              loadct,8,/silent
           end
        endcase     
        color = 249-200*wres/(max(wres))
        for i=0,xysize-1 do begin
           for j=0,xysize-1 do begin
              xfill = [a[i],a[i],a[i+1],a[i+1]]
              yfill = [a[j],a[j+1],a[j+1],a[j]]
              if(wres(i,j) gt 0) then begin
                 polyfill,xfill,yfill,color=color(i,j)
              endif   
           endfor
        endfor
     endfor
  endelse
end
