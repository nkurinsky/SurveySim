;Makes a histogram of the colors from model sources, saves the
;color arrays as an idl save file named model_hist
;Current as of 6/19

pro model_color,min,max,bin

  !p.thick=5
  !x.thick=5
  !y.thick=5
  !p.charthick=5
  !p.charsize=1.5

  set_plot,'ps'

  model = mrdfits('model.fits',1,/silent)
  bands = model.wavelength
  fluxes = model.models

  redshift = get_z_dist()
  mod_num = randomn(seed,n_elements(redshift))
  mod_num -= min(mod_num)
  mod_num *= 10.999/max(mod_num)
  
  t_num = mod_num*7.999/max(mod_num)
  
  mod_num = fix(mod_num)
  t_num = fix(t_num)
  
  num = n_elements(redshift)
  e250 = spire_errors(num,250)
  e350 = spire_errors(num,350)
  e500 = spire_errors(num,500)

  device,filename='plots/SPIRE_redshift_distribution.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  h = histogram(redshift,nbins=200,locations=xz)
  plot,xz,h,psym=10,xtitle=textoidl('z'),ytitle=textoidl('Number')
  device,/close

  device,filename='plots/SPIRE_model_distribution.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  h = histogram(mod_num,binsize=1,locations=xh)
  plot,xh,h,psym=10,xtitle=textoidl('\Beta Index'),ytitle=textoidl('Number')
  device,/close

  device,filename='plots/SPIRE_temperature_distribution.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  h = histogram(t_num,binsize=1,locations=xh)
  plot,xh,h,psym=10,xtitle=textoidl('Temperature Index'),ytitle=textoidl('Number')
  device,/close

  b250e = 250.d-6
  b350e = 350.d-6
  b500e = 500.d-6
  
  b1 = b250e/(1+redshift)
  b2 = b350e/(1+redshift)
  b3 = b500e/(1+redshift)
  f1 = fltarr(n_elements(redshift))
  f2 = fltarr(n_elements(redshift))
  f3 = fltarr(n_elements(redshift))

  for i=0,n_elements(redshift)-1 do begin
     bi = mod_num[i]
     ti = t_num[i]

     f1[i] = interpol(reform(fluxes(bi,ti,*)),bands,b1[i])+e250[i]
     f2[i] = interpol(reform(fluxes(bi,ti,*)),bands,b2[i])+e350[i]
     f3[i] = interpol(reform(fluxes(bi,ti,*)),bands,b3[i])+e500[i]
     ;;if((f1[i] le 0.d0) or (f2[i] le 0.d0) or (f3[i] le 0.d0))  then print,bi,ti,redshift[i]
  endfor

  gpts = where((f1 gt 0.d0) and (f2 gt 0.d0) and (f3 gt 0.d0))
  f1 = f1[gpts]
  f2 = f2[gpts]
  f3 = f3[gpts]
  
  c1 = acolor(f1,f3,b1,b3)
  c2 = acolor(f2,f3,b2,b3)

  res = histogram(c1,min=min,max=max,binsize=bin,locations=xc1)
  model_c1 = {hist:res,x:xc1}

  device,filename='plots/SPIRE_model_c1_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  plot,xc1,res/total(res),psym=10
  device,/close

  device,filename='plots/SPIRE_model_color_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color

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

  modelhist = {hist:res,num:n_elements(c1)}
  save,modelhist,model_c1,filename='save_files/model_hist.save'

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
