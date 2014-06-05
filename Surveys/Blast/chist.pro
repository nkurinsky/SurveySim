;Function for creating two dimensional histogram
;Current as of 6/22

function chist,x,y,xmax=xmax,ymax=ymax,xmin=xmin,ymin=ymin,binsize=binsize,xlocs=xlocs,ylocs=ylocs

  if((n_elements(x) eq 0) or (n_elements(y) eq 0)) then begin
     print,"Correct Calling Sequence is result = chist(x,y)"
     return,-1
  endif

  if(n_elements(xmax) eq 0) then xmax = max(x)
  if(n_elements(xmin) eq 0) then xmin = min(x)
  if(n_elements(ymax) eq 0) then ymax = max(y)
  if(n_elements(ymin) eq 0) then ymin = min(y)
  if(n_elements(binsize) eq 0) then binsize = 1

  xsize = ((xmax-xmin)/binsize)+1
  ysize = ((ymax-ymin)/binsize)+1

  xh = (findgen(xsize)/(xsize-1))*(xmax-xmin)+xmin
  yh = (findgen(ysize)/(ysize-1))*(ymax-ymin)+ymin

  res = lonarr(xsize-1,ysize-1)*0
  xlocs = fltarr(xsize-1)
  ylocs = fltarr(ysize-1)

  for xi=0,xsize-2 do begin
     gpts = where((x ge xh[xi]) and (x lt xh[xi+1]))
     if(gpts[0] ne -1) then begin
        ty = y[gpts]
        xlocs[xi] = xh[xi]+0.5*binsize
        for yi=0,ysize-2 do begin
           gpts = where((ty ge yh[yi]) and (ty lt yh[yi+1]))
           if (gpts[0] ge 0) then res(xi,yi) = n_elements(gpts)
           ylocs[yi] = yh[yi]+0.5*binsize
        endfor
     endif else begin
        xlocs[xi] = xh[xi]+0.5*binsize
        for yi=0,ysize-2 do begin
           res(xi,yi) = 0
           ylocs[yi] = yh[yi]+0.5*binsize
        endfor
     endelse
  endfor

  return,res

end
