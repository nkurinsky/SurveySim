pro simulation_results,savefile=savefile

  COMMON simulation_com

  if(n_elements(parameters) eq 0) then begin
     if(keyword_set(savefile)) then begin
        temp = load_parameters(savefile)
     endif else begin
        temp = load_parameters()
     endelse
     parameters = temp
  endif

  size_screen=get_screen_size()
  size_screen_alt = size_screen*0.85
  size_screen = size_screen*0.8
  gmain = widget_base(title='Simulation Output',/column,xsize=size_screen[0],ysize=size_screen_alt[1])
  r1 = widget_base(gmain,/row)
  r2 = widget_base(gmain,/row)
  r2b = widget_base(gmain,/row)
  r3 = widget_base(gmain,/row)

  widget_control,gmain,set_uvalue=mnum
  xdim = fix(size_screen[0]/3.0)
  ydim = fix(size_screen[1]/3.0)

  loadct,0,/silent

  lumfunct = widget_draw(r1,xsize=xdim,ysize=ydim)
  redshift = widget_draw(r1,xsize=xdim,ysize=ydim)
  models = widget_draw(r1,xsize=xdim,ysize=ydim)
  dcount1 = widget_draw(r2,xsize=xdim,ysize=ydim)
  dcount2 = widget_draw(r2,xsize=xdim,ysize=ydim)
  dcount3 = widget_draw(r2,xsize=xdim,ysize=ydim)
  sim_colors = widget_draw(r2b,xsize=xdim,ysize=ydim)
  obs_colors = widget_draw(r2b,xsize=xdim,ysize=ydim)
  comp_colors = widget_draw(r2b,xsize=xdim,ysize=ydim)
  
  diags = widget_button(r3,uvalue='diags',value='Diagnostics')
  refresh = widget_button(r3,uvalue='refresh',value='Refresh')
  close = widget_button(r3,uvalue='close',value='Close')
  quit = widget_button(r3,uvalue='quit',value='Quit')
  
  widget_control,gmain,/realize
  xmanager,'simulation_results',gmain,/no_block
  
  set_plot,'x'
  device,decomposed=0
  plot_settings,plot_type='x'
  loadct,0,/silent

  res = mrdfits(parameters.files.oname,0,head,/silent)
  dists = mrdfits(files.oname,3,head,/silent)

  gpts = where(dists.f3 gt 0)
  f1 = dists[gpts].f1
  f2 = dists[gpts].f2
  f3 = dists[gpts].f3
  z = dists[gpts].z
  lum = dists[gpts].lum

  xd1 = dists.s1
  xd2 = dists.s2
  xd3 = dists.s3
  yd1 = dists.dnds1
  yd2 = dists.dnds2
  yd3 = dists.dnds3

  gpts = where(yd1 gt 0)
  xd1 = xd1[gpts]/1.d3
  yd1 = yd1[gpts]

  gpts = where(yd2 gt 0)
  xd2 = xd2[gpts]/1.d3
  yd2 = yd2[gpts]

  gpts = where(yd3 gt 0)
  xd3 = xd3[gpts]/1.d3
  yd3 = yd3[gpts]

  counts = compute_counts(files.oname)

  pnum_out = fxpar(head,'tfields')

  widget_control,lumfunct,get_value=index
  wset,index
  device,decomposed=0
  plot_lumfunct

  widget_control,redshift,get_value=index
  wset,index
  h = histogram(z,binsize=0.1,locations=xh,min=0.2,max=5.0)
  plot,xh,h,psym=10,xrange=[0,max(z)],xstyle=1,xtitle='z',ytitle='dN/dz',title='Redshift Distribution'

  widget_control,models,get_value=index
  wset,index
  f = alog10(lum)
  h = histogram(f,locations=xh,min=min(f),max=max(f))
  xpts=10.0^(xh)
  plot,xpts,h,psym=10,/xlog,/ylog,xstyle=1,ystyle=0,title="Luminosity Distribution",xtitle="L",ytitle="dN/dL"

  widget_control,dcount1,get_value=index
  wset,index
  
  ;Herschel ATLAS counts at 250,350 and 500 (Clements et al. 2010)
  if( file_test('counts_clements10.dat')) then begin
     readcol,'counts_clements10.dat',skipline=2,numline=16,flux,nbin,corr,int_counts,int_err,diff_counts,diff_err,/silent
     flux /= 1.d3
     xrange=[min([xd1,flux]),max([xd1,flux])]
     yrange=[min([yd1,diff_counts,c1minus])/1.2,max([yd1,diff_counts,c1plus])*1.2]
     plot,flux,diff_counts,psym=1,symsize=2,xtitle=TeXtoIDL('F_{250}[Jy]'),ytitle=TeXtoIDL('(dN/dS)S^{2.5} [gal ster^{-1} J^{1.5}]'),/xlog,/ylog,title='Band 1 Counts',yrange=yrange,xrange=xrange,ystyle=1,xstyle=1
     oploterr,flux,diff_counts,diff_err
     oplot,xd1,yd1,psym=2
  endif else begin
     print,'Error: File "counts_clements10.dat" not found'
     plot,xd1,yd1,psym=2,symsize=2,xtitle=TeXtoIDL('F_{250}[Jy]'),ytitle=TeXtoIDL('(dN/dS)S^{2.5} [gal ster^{-1} J^{1.5}]'),/xlog,/ylog,title='Band 1 Counts'
  endelse
  oplot,xd1,c1plus,linestyle=1
  oplot,xd1,c1minus,linestyle=1

  widget_control,dcount2,get_value=index
  wset,index

  if( file_test('counts_clements10.dat')) then begin
     readcol,'counts_clements10.dat',skipline=19,numline=13,flux,nbin,corr,int_counts,int_err,diff_counts,diff_err,/silent
     flux /= 1.d3
     xrange=[min([xd2,flux]),max([xd2,flux])]
     yrange=[min([yd2,diff_counts,c2minus])/1.2,max([yd2,diff_counts,c2plus])*1.2]
     plot,flux,diff_counts,psym=1,symsize=2,xtitle=TeXtoIDL('F_{250}[Jy]'),ytitle=TeXtoIDL('(dN/dS)S^{2.5} [gal ster^{-1} J^{1.5}]'),/xlog,/ylog,title='Band 2 Counts',yrange=yrange,xrange=xrange,ystyle=1,xstyle=1
     oploterr,flux,diff_counts,diff_err
     oplot,xd2,yd2,psym=2
  endif else begin
     plot,xd2,yd2,psym=2,symsize=2,xtitle=TeXtoIDL('F_{250}[Jy]'),ytitle=TeXtoIDL('(dN/dS)S^{2.5} [gal ster^{-1} J^{1.5}]'),/xlog,/ylog,title='Band 2 Counts'
  endelse
  oplot,xd2,c2plus,linestyle=1
  oplot,xd2,c2minus,linestyle=1

  widget_control,dcount3,get_value=index
  wset,index

  if( file_test('counts_clements10.dat')) then begin
     readcol,'counts_clements10.dat',skipline=33,numline=10,flux,nbin,corr,int_counts,int_err,diff_counts,diff_err,/silent
     flux /= 1.d3
     xrange=[min([xd3,flux]),max([xd3,flux])]
     yrange=[min([yd3,diff_counts,c3minus])/1.2,max([yd3,diff_counts,c3plus])*1.2]
     plot,flux,diff_counts,psym=1,symsize=2,xtitle=TeXtoIDL('F_{250}[Jy]'),ytitle=TeXtoIDL('(dN/dS)S^{2.5} [gal ster^{-1} J^{1.5}]'),/xlog,/ylog,title='Band 3 Counts',yrange=yrange,xrange=xrange,ystyle=1,xstyle=1
     oploterr,flux,diff_counts,diff_err
     oplot,xd3,yd3,psym=2
  endif else begin
     plot,xd3,yd3,psym=2,symsize=2,xtitle=TeXtoIDL('F_{250}[Jy]'),ytitle=TeXtoIDL('(dN/dS)S^{2.5} [gal ster^{-1} J^{1.5}]'),/xlog,/ylog,title='Band 3 Counts'
  endelse
  oplot,xd3,c3plus,linestyle=1
  oplot,xd3,c3minus,linestyle=1

  comp = mrdfits(files.oname,0,head,/silent)
  model = mrdfits(files.oname,1,/silent)
  obs = mrdfits(files.oname,2,/silent)

  xysize = fxpar(head,'DIM')
  hist_min = fxpar(head,'H_MIN')
  hist_max = fxpar(head,'H_MAX')
  binsize = fxpar(head,'BINSIZE')
  
  a = findgen(xysize+1)
  a *= binsize
  a += hist_min
  
  loadct,39,/silent
  device,decomposed=0
  widget_control,sim_colors,get_value=index
  wset,index
  
  plot,[hist_min,hist_max],[hist_min,hist_max],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}'),title='Model Color Distribution'
  hb = 0.5

  color = 256*model/(max(model)*1.2)+30
  color[where(model eq 0)] = 0

  for i=0,xysize-1 do begin
     for j=0,xysize-1 do begin
        
        xfill = [a[i],a[i],a[i+1],a[i+1]]
        yfill = [a[j],a[j+1],a[j+1],a[j]]

        if(model(i,j) gt 0) then begin
           polyfill,xfill,yfill,color=color(i,j)
        endif

                                ;if(i eq xysize-1) then
                                ;oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
     endfor
     ;oplot,[a[i+1],a[i+1]],[hist_min,hist_max],linestyle=1
  endfor

  widget_control,obs_colors,get_value=index
  wset,index
  
  plot,[hist_min,hist_max],[hist_min,hist_max],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}'),title='Observed Color Distribution'
  hb = 0.5

  color = 256*obs/(max(obs)*1.2)+30
  color[where(obs eq 0)] = 0

  for i=0,xysize-1 do begin
     for j=0,xysize-1 do begin
        
        xfill = [a[i],a[i],a[i+1],a[i+1]]
        yfill = [a[j],a[j+1],a[j+1],a[j]]

        if(obs(i,j) gt 0) then begin
           polyfill,xfill,yfill,color=color(i,j)
        endif

                                ;if(i eq xysize-1) then
                                ;oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
     endfor
     ;oplot,[a[i+1],a[i+1]],[hist_min,hist_max],linestyle=1
  endfor

  widget_control,comp_colors,get_value=index
  wset,index
    
  pos = where(comp gt 0)
  neg = where(comp lt 0)
  resmod = comp
  resobs = abs(comp)
  if(neg[0] ne -1) then resmod[neg] = 0
  if(pos[0] ne -1) then resobs[pos] = 0

  plot,[hist_min,hist_max],[hist_min,hist_max],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}'),title='Color Distribution Comparison'
  hb = 0.5

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
           
                                ;if(i eq xysize-1) then
                                ;oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
        endfor
        ;oplot,[a[i+1],a[i+1]],[hist_min,hist_max],linestyle=1
     endfor
  endfor

end
