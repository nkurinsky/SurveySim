;Makes a histogram of the colors from observed sources, saves the
;color arrays as an idl save file named obs_hist
;Current as of 6/22

pro obs_color,min,max,bin
  
  !p.thick=5
  !x.thick=5
  !y.thick=5
  !p.charthick=5
  !p.charsize=1.5

  restore,'save_files/matched_fluxes.save'
  c1 = acolor(f250,f500,250.d0,500.d0)
  c2 = acolor(f350,f500,350.d0,500.d0)

  res = histogram(c1,min=min,max=max,binsize=bin,locations=xc1)
  obs_c1 = {hist:res,x:xc1}

  device,filename='plots/obs_c1_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  plot,xc1,res/total(res),psym=10
  device,/close

  loadct,1,/silent

  set_plot,'ps'
  device,filename='plots/obs_color.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman'
  plot,c1,c2,psym=2,xrange=[min,max],yrange=[min,max]

  device,/close

  device,filename='plots/obs_color_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color

  max1 = max
  max2 = max
  min1 = min
  min2 = min
  binsize = bin

  res = chist(c1,c2,xmax=max1,ymax=max2,xmin=min1,ymin=min2,binsize=binsize,xlocs=x,ylocs=y)
  
  plot,[min1,max1],[min2,max2],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}')
  hb = binsize*0.5
  scale = ceil(299/max(res))

  color = 260*res/(max(res)*1.2)+30
  color[where(res eq 0)] = 0

  obshist = {hist:res,num:n_elements(c1)}
  save,obshist,obs_c1,filename='save_files/obs_hist.save'

  loadct,39,/silent

  for i=0,n_elements(x)-1 do begin
     
     for j=0,n_elements(y)-1 do begin
        
        xfill = [x[i]-hb,x[i]-hb,x[i]+hb,x[i]+hb]
        yfill = [y[j]-hb,y[j]+hb,y[j]+hb,y[j]-hb]

        if(res(i,j) gt 0) then begin
           polyfill,xfill,yfill,color=color(i,j)
        endif

        if(i eq n_elements(x)-1) then oplot,[min1,max2],[y[j]+0.5*binsize,y[j]+0.5*binsize],linestyle=1        
     endfor
     oplot,[x[i]+0.5*binsize,x[i]+0.5*binsize],[min2,max2],linestyle=1
  endfor

  device,/close

end
