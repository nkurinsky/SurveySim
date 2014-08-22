pro read_output,filename
  
  print,filename
  dists = mrdfits(filename,3,head,/silent)

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

  count_dists = mrdfits(filename,6,head,/silent)

  alpha = 0.159                ;one std deviation
  plusfrac = (1.0-alpha)
  minusfrac = alpha

  chis = count_dists.chisq
  gpts = where(chis gt median(chis))
  count_dists = count_dists[gpts]
  c1=count_dists.dnds250
  c2=count_dists.dnds350
  c3=count_dists.dnds500

  c1size = n_elements(count_dists[0].dnds250)
  c2size = n_elements(count_dists[0].dnds350)
  c3size = n_elements(count_dists[0].dnds500)

  c1mean = make_array(c1size,value=0.0)
  c2mean = make_array(c2size,value=0.0)
  c3mean = make_array(c3size,value=0.0)
  c1plus = make_array(c1size,value=0.0)
  c2plus = make_array(c2size,value=0.0)
  c3plus = make_array(c3size,value=0.0)
  c1minus = make_array(c1size,value=0.0)
  c2minus = make_array(c2size,value=0.0)
  c3minus = make_array(c3size,value=0.0)

    for i=0,c1size-1 do begin
     dnds = c1[c1size-i-1,*]
     dnds = dnds[sort(dnds)]
     dnds = dnds[where(dnds gt 0)]
     pi = plusfrac*n_elements(dnds)
     mi = minusfrac*n_elements(dnds)
     c1mean[i] = mean(dnds)
     c1plus[i] = dnds[pi]
     c1minus[i] = dnds[mi]
  endfor

  for i=0,c2size-1 do begin
     dnds = c2[i,*]
     dnds = dnds[sort(dnds)]
     dnds = dnds[where(dnds gt 0)]
     pi = plusfrac*n_elements(dnds)
     mi = minusfrac*n_elements(dnds)
     c2mean[i] = mean(dnds)
     c2plus[i] = dnds[pi]
     c2minus[i] = dnds[mi]
  endfor

  for i=0,c3size-1 do begin
     dnds = c3[i,*]
     dnds = dnds[sort(dnds)]
     dnds = dnds[where(dnds gt 0)]
     pi = plusfrac*n_elements(dnds)
     mi = minusfrac*n_elements(dnds)
     c3mean[i] = mean(dnds)
     c3plus[i] = dnds[pi]
     c3minus[i] = dnds[mi]
  endfor

  set_plot,'ps'
  plot_settings,plot_type='ps'
  device,filename='redshift_dist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  h = histogram(z,binsize=0.1,locations=xh,min=0.2,max=5.0)
  plot,xh,h,psym=10,xrange=[0,max(z)],xstyle=1,xtitle='z',ytitle='dN/dz',title='Redshift Distribution'
  device,/close

  device,filename='band1_counts.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
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
  device,/close
  
  device,filename='band2_counts.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
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
  device,/close

  device,filename='band3_counts.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
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
  device,/close

  comp = mrdfits(filename,0,head,/silent)
  model = mrdfits(filename,1,/silent)
  obs = mrdfits(filename,2,/silent)

  xysize = fxpar(head,'DIM')
  hist_min = fxpar(head,'H_MIN')
  hist_max = fxpar(head,'H_MAX')
  binsize = fxpar(head,'BINSIZE')
  
  alpha = fxpar(head,'ALPHA')
  beta = fxpar(head,'BETA')
  phi0 = fxpar(head,'PHI0')
  L0 = fxpar(head,'L0')
  p = fxpar(head,'P')
  q = fxpar(head,'Q')
  zcut = fxpar(head,'ZCUT')

  a = findgen(xysize+1)
  a *= binsize
  a += hist_min

  device,filename='model_color_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  
  plot,[hist_min,hist_max],[hist_min,hist_max],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}')
  hb = 0.5

  color = 260*model/(max(model)*1.2)+30
  color[where(model eq 0)] = 0

  loadct,39,/silent

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

  device,/close

  device,filename='obs_color_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  
  plot,[hist_min,hist_max],[hist_min,hist_max],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}')
  hb = 0.5

  color = 260*obs/(max(obs)*1.2)+30
  color[where(obs eq 0)] = 0

  loadct,39,/silent

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

  device,/close

  device,filename='comp_color_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  
  pos = where(comp gt 0)
  neg = where(comp lt 0)
  resmod = comp
  resobs = abs(comp)
  if(neg[0] gt -1) then resmod[neg] = 0
  if(pos[0] gt -1) then resobs[pos] = 0

  plot,[hist_min,hist_max],[hist_min,hist_max],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}')
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

  device,/close

  omega_m = 0.28
  omega_l = 0.72
  
  device,filename='sim_lumfunct.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated

  plot,[8,13],[1e-10,1e0],/ylog,/nodata,xtitle=textoidl('log_{10}(L_{fir}) [L_{sun}]'),ytitle=textoidl('log_{10}(N(L_{fir})/\Omega dV_c)'),ystyle=1
  lums = indgen(21)/4.0+8.0
  
  dz = 0.25
  loadct,39,/silent
  for z=dz,5.0,dz do begin
     if z le zcut then begin
        t1 = (10^phi0)*((1+z)^p)
        t2 = (10^L0)*((1+z)^q)
     endif else begin
        t1 = (10^phi0)*((3.0)^p)
        t2 = (10^L0)*((3.0)^q)
     endelse
     r = 10^lums/t2
     nsrcs=t1/(r^alpha+r^beta)
     oplot,lums,nsrcs,color=z*40+20
  endfor

  for i=0.0,5.0,0.1 do begin
     ind = (9+2*i/5)
     oplot,[ind,ind+0.02,ind,ind+0.02,ind,ind+0.02],[10^(-8),10^(-8),10^(-7.9),10^(-7.9),10^(-7.8),10^(-7.8)],psym=2,color=i*40+20,symsize=0.75
  endfor

  loadct,0,/silent
  oplot,[8,13],[1,1],linestyle=1

  xyouts,9.0,10^(-7.2),textoidl('Redshift')
  xyouts,9.0,10^(-7.6),textoidl('0')
  xyouts,10.9,10^(-7.6),textoidl('5')

  device,/close

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
  accept = make_array(n_elements(apts),n_elements(res.(cpts[0])),value=0.0)
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

  set_plot,'ps'
  plot_settings,plot_type='ps'
  device,filename='fit_results.eps',xsize=10,ysize=8,/inches,/times,set_font='Times-Roman',/color,/encapsulated
  multiplot,[dim,dim],/init,/rowmajor,mTitle="MCMC Fitting Results",gap=0.005
  cgText, 0.6, 0.9, Alignment=0, /Normal, 'Fitting Results:', Charsize=1.25

  for x=0,dim-1 do begin
     
     ppts = where(strmatch(alltags,tags[x]+'*',/FOLD_CASE) eq 1)
     p = make_array(n_elements(ppts),n_elements(res.(ppts[0])),value = 0.0)
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
           cgText, 0.6, 0.87-0.03*x, Alignment=0, /Normal, Charsize=1.25,textoidl(tags[x]+": "+strcompress(string(stats[0],format='(D0.3)'))+"\pm"+strcompress(string(stats[1],format='(D0.3)')))
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
                    polyfill,xfill,yfill,color=color_scale*(alog(chi_plot))
                 endif
              endfor
           endfor

        endif 
     endfor
  endfor
  
  device,/close
  
  multiplot,/reset

  gpts = where(accept eq 1.0)
  crange=[min(chis[gpts]),max(chis[gpts])]
  n = n_elements(res)

  device,filename='chisq_v_run.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated

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
  
  device,/close

  device,filename='chisq_hist.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated
  
  gpts = where(accept eq 1.0)
  histrange = [0,min(chis)*10]
  chist = histogram(chis,nbins=50,locations=xchist,min=histrange[0],max=histrange[1])
  chist_acpt = histogram(chis[gpts],nbins=50,locations=xchist_acpt,min=histrange[0],max=histrange[1])
  cmax = max(chist)
  gpts = where(chist lt 1)
  chist[gpts] = 0.001
  plot,xchist,chist,psym=10,xstyle=1,ystyle=1,xrange=histrange,yrange=[0.9,cmax^1.2],xtitle=Textoidl("\chi^2"),ytitle="N",title=textoidl("Total \chi^2 Distribution"),/ylog
  oplot,xchist_acpt,chist_acpt,psym=10,linestyle=1

  device,/close

  device,filename='convergence.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated

  res= mrdfits(filename,5,/silent)
  gpts = where(res.(0) gt 0)
  res = res[gpts]
  rmax = 0
  for i=0,n_elements(tag_names(res))-1 do begin
     tmax = max(res.(i))
     if (tmax gt rmax) then rmax = tmax
  endfor
  
  plot,[0,n_elements(res)],[1.0,rmax],xstyle=1,ystyle=1,title="Convergence",xtitle="Test Number",ytitle="R",/nodata
  ;taken out of the program, this is invalid...maybe have optional model fits file
  ;oplot,[0,n_elements(gpts)],[msettings.conv_rmax,msettings.conv_rmax]
  for i=0,n_elements(tag_names(res))-1 do begin
     oplot,res.(i),linestyle=i
  endfor
  
  device,/close

  if(file_test('plots')) then begin
     file_delete,'plots',/recursive
  endif
  file_mkdir,'plots'
  file_move,'*.eps','plots'

end
