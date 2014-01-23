
;==================================================================
; Writen by Noah Kurinsky, version recent as of 10/17/13
; Edited by Anna Sajina, December 2012
; certain files in the same directory are required for proper
; function of this widget:
;
;==================================================================

pro simulation_output

  COMMON simulation_com,info,ldata,ldat0,sdat,cdat,bands,msettings,files
  
  !p.thick=5
  !x.thick=5
  !y.thick=5
  !p.charthick=5
  !p.charsize=1.5
  
  comp = mrdfits(files.oname,0,head,/silent)
  model = mrdfits(files.oname,1,/silent)
  obs = mrdfits(files.oname,2,/silent)

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

  set_plot,'ps'
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

        if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
     endfor
     oplot,[a[i+1],a[i+1]],[hist_min,hist_max],linestyle=1
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

        if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
     endfor
     oplot,[a[i+1],a[i+1]],[hist_min,hist_max],linestyle=1
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
           
           if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
        endfor
        oplot,[a[i+1],a[i+1]],[hist_min,hist_max],linestyle=1
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

  res = mrdfits(files.oname,0,head,/silent)
  alpha = fxpar(head,'ALPHA')
  beta = fxpar(head,'BETA')
  phi0 = fxpar(head,'PHI0')
  L0 = fxpar(head,'L0')
  p = fxpar(head,'P')
  q = fxpar(head,'Q')

  loadct,0,/silent
  
; chain read operations 
  res = mrdfits(files.oname,4,/silent)
  tags = tag_names(res)
  p = [res.p0,res.p1,res.p2,res.p3,res.p4]
  q = [res.q0,res.q1,res.q2,res.q3,res.q4]
  chis = [res.chisq0,res.chisq1,res.chisq2,res.chisq3,res.chisq4]
  prange=[min(p),max(p)]
  qrange=[min(q),max(q)]
  crange=[min(chis)/1.2,min(chis)*10]
  log_crange=alog(crange)
  dp=(prange[1]-prange[0])/50.0
  dq=(qrange[1]-qrange[0])/50.0
  dc=alog(crange[1]-crange[0])/50.0
  n = n_elements(res)

  device,filename='pspace.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated

  plot,prange,qrange,xstyle=1,ystyle=1,/nodata,xtitle='P',ytitle='Q',title="Chain Coverage of Parameter Space"
  loadct,39,/silent
  p1 = res.p0
  p2 = res.q0
  oplot,p1,p2,psym=2,symsize=0.25,color=40
  xyouts,p1[0]+0.1,p2[0]+0.1,"Start"
  xyouts,p1[n-1]+0.1,p2[n-1]+0.1,"End"
  p1 = res.p1
  p2 = res.q1
  oplot,p1,p2,psym=2,symsize=0.25,color=80
  xyouts,p1[0]+0.1,p2[0]+0.1,"Start"
  xyouts,p1[n-1]+0.1,p2[n-1]+0.1,"End"
  p1 = res.p2
  p2 = res.q2
  oplot,p1,p2,psym=2,symsize=0.25,color=120
  xyouts,p1[0]+0.1,p2[0]+0.1,"Start"
  xyouts,p1[n-1]+0.1,p2[n-1]+0.1,"End"
  p1 = res.p3
  p2 = res.q3
  oplot,p1,p2,psym=2,symsize=0.25,color=160
  xyouts,p1[0]+0.1,p2[0]+0.1,"Start"
  xyouts,p1[n-1]+0.1,p2[n-1]+0.1,"End"
  p1 = res.p4
  p2 = res.q4
  oplot,p1,p2,psym=2,symsize=0.25,color=200
  xyouts,p1[0]+0.1,p2[0]+0.1,"Start"
  xyouts,p1[n-1]+0.1,p2[n-1]+0.1,"End"

  device,/close
  device,filename='chisq_v_run.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated

  loadct,1,/silent
  plot,[0,n_elements(p1)],crange,/ylog,xstyle=1,ystyle=1,/nodata,xtitle='Run Number',ytitle=textoidl("\chi^2"),title="Temporal Likelihood Trends"
  loadct,39,/silent
  oplot,res.chisq0,color=40
  oplot,res.chisq1,color=80
  oplot,res.chisq2,color=120
  oplot,res.chisq3,color=160
  oplot,res.chisq4,color=200
  
  device,/close
  device,filename='chisq_v_p.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated

  loadct,1,/silent
  plot,prange,crange,/ylog,xstyle=1,ystyle=1,/nodata,xtitle='P',ytitle=textoidl("\chi^2"),title=textoidl("P \chi^2 Distribution")
  loadct,39,/silent
  
  oplot,res.p0,res.chisq0,color=40,psym=2,symsize=0.25
  oplot,res.p1,res.chisq1,color=80,psym=2,symsize=0.25
  oplot,res.p2,res.chisq2,color=120,psym=2,symsize=0.25
  oplot,res.p3,res.chisq3,color=160,psym=2,symsize=0.25
  oplot,res.p4,res.chisq4,color=200,psym=2,symsize=0.25
  
  device,/close

  device,filename='mean_likelihood.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated

  plot,prange,qrange,/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle="P",ytitle="Q",title='Mean Likelihood Space'
  chi_min = alog(min(chis))
  chi_max = alog(max(chis))
  color_scale = 256/((chi_max-chi_min)*1.2)
  
  for i=prange[0],prange[1]-dp/2,dp do begin
     for j=qrange[0],qrange[1]-dq/2,dq do begin
        
        xfill = [i,i,i+dp,i+dp]
        yfill = [j,j+dq,j+dq,j]

        gpts = where((p gt i) and (p le i+dp) and (q gt j) and (q lt j+dq))
        chi_mean = mean(chis[gpts])

        if(gpts[0] ne -1) then begin
           polyfill,xfill,yfill,color=color_scale*(chi_max - alog(chi_mean))
        endif
        
        ;if(i ge prange[1]-dp) then oplot,prange,[j+dq,j+dq],linestyle=1        
     endfor
     ;oplot,[i+dp,i+dp],qrange,linestyle=1
  endfor

  device,/close

  device,filename='convergence.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated
  
  res= mrdfits(files.oname,5,/silent)
  gpts = where(res.(0) gt 0)
  res = res[gpts]
  rmax = 0
  for i,n_elements(tag_names(res))-1 do begin
     tmax = max(res.(rpts[i]))
     if (tmax gt rmax) then rmax = tmax
  endfor
  
  plot,[0,n_elements(res)],[1.0,rmax],xstyle=1,ystyle=1,title="Convergence",xtitle="Test Number",ytitle="R",/nodata
  oplot,[0,n_elements(gpts)],[msettings.conv_rmax,msettings.conv_rmax]
  for i,n_elements(rpts)-1 do begin
     oplot,res.(i),linestyle=i
  endfor

  device,/close

  if(file_test('plots')) then begin
     file_delete,'plots',/recursive
  endif
  file_mkdir,'plots'
  file_move,'*.eps','plots'

end
