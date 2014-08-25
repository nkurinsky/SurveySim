pro simulation_diagnostics

  COMMON simulation_com

  size_screen=get_screen_size()
  size_screen_alt = size_screen*0.85
  size_screen = size_screen*0.8
  dmain = widget_base(title='Simulation Diagnostics',/column,xsize=size_screen[0],ysize=size_screen_alt[1])
  r1 = widget_base(dmain,/row)
  c1 = widget_base(r1,/column)
  r3 = widget_base(dmain,/row)

  widget_control,dmain,set_uvalue=mnum
  xdim = fix(size_screen[0]/3.0)
  ydim = fix(size_screen[1]/3.0)

  res = mrdfits(files.oname,0,head,/silent)
  alpha = fxpar(head,'ALPHA')
  beta = fxpar(head,'BETA')
  phi0 = fxpar(head,'PHI0')
  L0 = fxpar(head,'L0')
  p = fxpar(head,'P')
  q = fxpar(head,'Q')

  loadct,0,/silent

  chisqr = widget_draw(c1,xsize=xdim,ysize=ydim)
  chidist = widget_draw(c1,xsize=xdim,ysize=ydim)
  conv = widget_draw(c1,xsize=xdim,ysize=ydim)
  resplot = widget_draw(r1,xsize=2*xdim,ysize=3*ydim)
  
  graphs = widget_button(r3,uvalue='graphs',value='Return to Output')
  refresh = widget_button(r3,uvalue='refresh',value='Refresh')
  close = widget_button(r3,uvalue='close',value='Close')
  quit = widget_button(r3,uvalue='quit',value='Quit')

  widget_control,dmain,/realize
  xmanager,'simulation_diagnostics',dmain,/no_block
  
  set_plot,'x'
  plot_settings,plot_type='x'
  device,decomposed=0
  loadct,0,/silent
  
; chain read operations 
  res = mrdfits('output.fits',4,/silent)
  alltags = tag_names(res)
  tags = alltags[where(strmatch(alltags,'*0',/FOLD_CASE) eq 1)]

  for i=0,n_elements(tags)-1 do begin
     tags[i] = strmid(tags[i],0,strlen(tags[i])-1)
  endfor
  tags = tags[where((tags ne 'CHISQ') and (tags ne 'ACPT'))]
  print,'Fitted Variables: ',tags

  dim = n_elements(tags)
  nbins = 50.0
  
  cpts = where(strmatch(alltags,'CHISQ*',/FOLD_CASE) eq 1)
  for ci=0,n_elements(cpts)-1 do begin
     if(ci eq 0) then begin
        chis = res.(cpts[ci])
     endif else begin
        chis = [chis,res.(cpts[ci])]
     endelse
  endfor

  apts = where(strmatch(alltags,'ACPT*',/FOLD_CASE) eq 1)
  for ai=0,n_elements(apts)-1 do begin
     if(ai eq 0) then begin
        accept = res.(apts[ai])
     endif else begin
        accept = [accept,res.(apts[ai])]
     endelse
  endfor
  
  chi_med = median(chis)
  print,"median: ",chi_med
  chistpts = where(chis lt chi_med)

  chi_min = alog(min(chis))
  chi_max = alog(max(chis))
  color_scale = 256/((chi_max-chi_min)*1.2)

  widget_control,resplot,get_value=index
  wset,index

  multiplot,[dim,dim],/init,/rowmajor,mTitle="MCMC Fitting Results",gap=0.005
  cgText, 0.6, 0.9, Alignment=0, /Normal, 'Fitting Results:', Charsize=1.25,color='black'

  for x=0,dim-1 do begin
     ppts = where(strmatch(alltags,tags[x]+'*',/FOLD_CASE) eq 1)
     for pi=0,n_elements(ppts)-1 do begin
        if(pi eq 0) then begin
           p = res.(ppts[pi])
        endif else begin
           p = [p,res.(ppts[pi])]
        endelse
     endfor
     prange = [min(p),max(p)]
     for y=0,dim-1 do begin
        multiplot
        if (x eq y) then begin  ;make histogram
           stats = moment(p[chistpts],maxmoment=2)
           print,"Variable: ",tags[x],", Mean: ",stats[0]," Variance: ",stats[1]
           cgText, 0.6, 0.87-0.03*x, Alignment=0, /Normal, Charsize=1.25,textoidl(tags[x]+": "+strcompress(string(stats[0],format='(D0.3)'))+"\pm"+strcompress(string(stats[1],format='(D0.3)'))),color='black'
           h = histogram(p[chistpts],locations=xh,max=prange[1],min=prange[0],nbins=nbins)
           h = float(h)/float(max(h))
           loadct,1,/silent
           
           xt = ''
           yt = ''
           if(x eq 0) then yt=tags[y]
           if(y eq dim-1) then xt=tags[x]
           plot,xh,h,psym=10,xrange=prange,yrange=[0,1.2],xtitle=xt,ytitle=yt,xstyle=1,ystyle=1

        endif else if (x lt y) then begin ;make liklihood space

           qpts = where(strmatch(alltags,tags[y]+'*',/FOLD_CASE) eq 1)
           for qi=0,n_elements(qpts)-1 do begin
              if(qi eq 0) then begin
                 q = res.(qpts[qi])
              endif else begin
                 q = [q,res.(qpts[qi])]
              endelse
           endfor
           qrange = [min(q),max(q)]
           dp=(prange[1]-prange[0])/nbins
           dq=(qrange[1]-qrange[0])/nbins

           xt = ''
           yt = ''
           if(x eq 0) then yt=tags[y]
           if(y eq dim-1) then xt=tags[x]
           
           loadct,1,/silent
           plot,prange,qrange,/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=xt,ytitle=yt
           loadct,39,/silent
           for i=prange[0],prange[1]-dp/2,dp do begin
              for j=qrange[0],qrange[1]-dq/2,dq do begin
                 
                 xfill = [i,i,i+dp,i+dp]
                 yfill = [j,j+dq,j+dq,j]
                 
                 gpts = where((p gt i) and (p le i+dp) and (q gt j) and (q lt j+dq))
                 chi_plot = n_elements(gpts)
                 
                 if(gpts[0] ne -1) then begin
                    polyfill,xfill,yfill,color=color_scale*alog(chi_plot)
                 endif
              endfor
           endfor

        endif 
     endfor
  endfor
  
  multiplot,/reset

  widget_control,chisqr,get_value=index
  wset,index

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
  
  widget_control,chidist,get_value=index
  wset,index

  gpts = where(accept eq 1.0)
  histrange = [0,min(chis)*10]
  chist = histogram(chis,nbins=50,locations=xchist,min=histrange[0],max=histrange[1])
  chist_acpt = histogram(chis[gpts],nbins=50,locations=xchist_acpt,min=histrange[0],max=histrange[1])
  cmax = max(chist)
  gpts = where(chist lt 1)
  chist[gpts] = 0.001
  plot,xchist,chist,psym=10,xstyle=1,ystyle=1,xrange=histrange,yrange=[0.9,cmax^1.2],xtitle=Textoidl("\chi^2"),ytitle="N",title=textoidl("Total \chi^2 Distribution"),/ylog
  oplot,xchist_acpt,chist_acpt,psym=10,linestyle=1

  widget_control,conv,get_value=index
  wset,index

  res= mrdfits(files.oname,5,/silent)
  gpts = where(res.(0) gt 0)
  res = res[gpts]
  rmax = 0
  for i=0,n_elements(tag_names(res))-1 do begin
     tmax = max(res.(i))
     if (tmax gt rmax) then rmax = tmax
  endfor
  
  plot,[0,n_elements(res)],[1.0,rmax],xstyle=1,ystyle=1,title="Convergence",xtitle="Test Number",ytitle="R",/nodata
  oplot,[0,n_elements(gpts)],[msettings.conv_rmax,msettings.conv_rmax]
  for i=0,n_elements(tag_names(res))-1 do begin
     oplot,res.(i),linestyle=i
  endfor

end
