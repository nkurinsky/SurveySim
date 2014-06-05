;Compares the colors from the model and observed save files as a
;histogram
;Current as of 6/19

pro compare,min,max,bin
  
  !p.thick=5
  !x.thick=5
  !y.thick=5
  !p.charthick=5
  !p.charsize=1.5

  restore,'save_files/model_hist.save'
  restore,'save_files/obs_hist.save'

  max1 = max
  max2 = max
  min1 = min
  min2 = min
  binsize = bin

  ;turn into percent of sample
  pmodel = (modelhist.hist/float(modelhist.num))
  pobs = (obshist.hist/float(obshist.num))
  res = pmodel-pobs
  pos = where(res gt 0)
  neg = where(res lt 0)
  resmod = res
  resobs = abs(res)
  resmod[neg] = 0
  resobs[pos] = 0

  ;c1 histogram
  c1hist_mod = model_c1.hist
  c1hist_mod /= total(c1hist_mod)
  c1hist_obs = obs_c1.hist
  c1hist_obs /= total(c1hist_obs)
  
  gpts = where((c1hist_mod gt 0.0) or (c1hist_obs gt 0.0))
  diffc1 = fltarr(n_elements(c1hist_mod))
  diffc1[gpts] = (c1hist_mod[gpts]-c1hist_obs[gpts])*100d0

  x = model_c1.x

  set_plot,'ps'
  crange = [x[0],x[n_elements(x)-1]]
  device,filename='plots/c1_comp_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  plot,x,diffc1,psym=10,yrange = [-25,25],xrange=crange,ystyle=1,xstyle=1
  oplot,crange,[0,0],linestyle=1
  device,/close

  ;dummy call for getting x and y
  temp = chist([0,0],[0,0],xmax=max1,ymax=max2,xmin=min1,ymin=min2,binsize=binsize,xlocs=x,ylocs=y)

  device,filename='plots/comparison_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color

  plot,[min1,max1],[min2,max2],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}')
  hb = binsize*0.5

  for resi=0,1 do begin

     case resi of
        0:begin
           wres = resmod
           loadct,1,/silent
        end
        1:begin
           wres = resobs
           loadct,3,/silent
        end
     endcase

     color = 249-200*wres/(max(abs(wres)))
     
     for i=0,n_elements(x)-1 do begin
        
        for j=0,n_elements(y)-1 do begin
           
           xfill = [x[i]-hb,x[i]-hb,x[i]+hb,x[i]+hb]
           yfill = [y[j]-hb,y[j]+hb,y[j]+hb,y[j]-hb]
           
           if(wres(i,j) gt 0) then begin
              polyfill,xfill,yfill,color=color(i,j)
           endif
           
           if((i eq n_elements(x)-1) and (resi eq 1)) then oplot,[min1,max2],[y[j]+0.5*binsize,y[j]+0.5*binsize],linestyle=1        
        endfor
        if (resi eq 1) then oplot,[x[i]+0.5*binsize,x[i]+0.5*binsize],[min2,max2],linestyle=1
     endfor
  endfor

  device,/close

end
