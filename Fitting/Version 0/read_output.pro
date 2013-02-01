pro read_output,do_command=do_command,command=command

  if(n_elements(command) eq 0) then command = 'test_fit'
  if(n_elements(do_command) eq 0) then do_command=0

  if(do_command) then spawn,command
  
  comp = mrdfits("output.fits",0,head)
  model = mrdfits("output.fits",1)
  obs = mrdfits("output.fits",2)
  dists = mrdfits("output.fits",3)

  xysize = fxpar(head,'DIM')

  set_plot,'ps'
  device,filename='model_color_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  
  plot,[0,xysize],[0,xysize],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1;,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}')
  hb = 0.5

  color = 260*model/(max(model)*1.2)+30
  color[where(model eq 0)] = 0

  loadct,39,/silent

  for i=0,xysize-1 do begin
     for j=0,xysize-1 do begin
        
        xfill = [i,i,i+1,i+1]
        yfill = [j,j+1,j+1,j]

        if(model(i,j) gt 0) then begin
           polyfill,xfill,yfill,color=color(i,j)
        endif

        if(i eq xysize-1) then oplot,[0,xysize],[j+1,j+1],linestyle=1        
     endfor
     oplot,[i+1,i+1],[0,xysize],linestyle=1
  endfor

  device,/close

  device,filename='obs_color_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  
  plot,[0,xysize],[0,xysize],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1;,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}')
  hb = 0.5

  color = 260*obs/(max(obs)*1.2)+30
  color[where(obs eq 0)] = 0

  loadct,39,/silent

  for i=0,xysize-1 do begin
     for j=0,xysize-1 do begin
        
        xfill = [i,i,i+1,i+1]
        yfill = [j,j+1,j+1,j]

        if(obs(i,j) gt 0) then begin
           polyfill,xfill,yfill,color=color(i,j)
        endif

        if(i eq xysize-1) then oplot,[0,xysize],[j+1,j+1],linestyle=1        
     endfor
     oplot,[i+1,i+1],[0,xysize],linestyle=1
  endfor

  device,/close

  device,filename='comp_color_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  
  plot,[0,xysize],[0,xysize],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1;,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}')
  hb = 0.5

  color = 260*comp/(max(comp)*1.2)+30
  color[where(comp eq 0)] = 0

  loadct,39,/silent

  for i=0,xysize-1 do begin
     for j=0,xysize-1 do begin
        
        xfill = [i,i,i+1,i+1]
        yfill = [j,j+1,j+1,j]

        if(comp(i,j) gt 0) then begin
           polyfill,xfill,yfill,color=color(i,j)
        endif

        if(i eq xysize-1) then oplot,[0,12],[j+1,j+1],linestyle=1        
     endfor
     oplot,[i+1,i+1],[0,xysize],linestyle=1
  endfor

  device,/close

  set_plot,'x'

  f1 = dists.f1
  f2 = dists.f2
  f3 = dists.f3
  z = dists.z
  m = dists.m

  gpts = where(f3 gt 0)
  f1 = f1[gpts]
  f2 = f2[gpts]
  f3 = f3[gpts]
  z = z[gpts]
  m = m[gpts]

  window,0
  h = histogram(z,nbins=50,locations=xh)
  plot,xh,h,psym=10,xrange=[0,10],xstyle=1

  window,1
  h = histogram(m,nbins=50,locations=xh)
  plot,xh,h,psym=10

  window,2
  h = histogram(f1,nbins=50,locations=xh)
  plot,xh,h,psym=10,/xlog,/ylog,yrange=[1e0,1e4],ystyle=1,xstyle=1
  
  window,3
  h = histogram(f2,nbins=50,locations=xh)
  plot,xh,h,psym=10,/xlog,/ylog,yrange=[1e0,1e4],ystyle=1,xstyle=1
  
  window,4
  h = histogram(f3,nbins=50,locations=xh)
  plot,xh,h,psym=10,/xlog,/ylog,yrange=[1e0,1e4],ystyle=1,xstyle=1
end
