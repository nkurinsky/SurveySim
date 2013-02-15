;==================================================================
; Writen by Noah Kurinsky, version recent as of 2/15/13
; Edited by Anna Sajina, December 2012
; certain files in the same directory are required for proper
; function of this widget:
;
;==================================================================

pro simulation

  COMMON simulation_com,cdir,info,ot,ldist,lum,settings,bands

;==================================================================
;the directory where the fitting code lives
;this should be generalized
;------------------------------------------------------------------
cdir='/Users/annie/students/noah_kurinsky/Fitting/v4/'
  
;==================================================================
;the INFO structure holds the key widget control parameters as well as 
;the basic simulation settings
;------------------------------------------------------------------
  info={$
;widget settings
       base: 0L, $
       base1: 0L, $
       base2: 0L, $
       base3:0L, $
       base_out: 0L, $
       base1_out: 0L, $
       base2_out: 0L, $
       base4:0L, $
       base5:0L, $
       xsize1:150., $
       ysize1:100., $
       xsize2:150., $
       ysize2:100., $
       magnification:1.0, $  
       draw1:0L, $
       draw2:0L, $
       draw3:0L, $
       draw4:0L, $
       draw5:0L, $
       draw6:0L, $     
       win_id1:0L, $
       win_id2:0L, $
       win_id3:0L, $
       ncolors: 0L, $
;basic simulation settings
       ofile:'observation.save', $
       mfile:'model.save', $
       oname:'output.fits', $
       pnum:2,$
       znum:2 }
;====================================================================
; Colors
;--------------------------------------------------------------------
;device, truecolor ;, depth=24 ;decompose=0
  loadct,3,/silent
  tvlct,red,green,blue,/get
  info.ncolors=!d.TABLE_SIZE
  red[info.ncolors-1]=0B
  green[info.ncolors-1]=255B
  blue[info.ncolors-1]=0B
  tvlct,red,green,blue
  
;; Screen size
  size_screen = get_screen_size()
  size_screen=size_screen*0.8
  spectrum_xsize=size_screen(0)*2./3. & spectrum_ysize=size_screen(1)/5.
  
  plotSize = 100.
  info.magnification = size_screen[1]/(plotSize+info.ysize1)*0.9
  
  ofile = info.ofile
  mfile = info.mfile
  pnum = long(info.pnum)
  znum = long(info.znum)

;;widget base initialization
  info.base = widget_base(title='Model Setup and Initialization',/column) ;main base
  info.base2= widget_base(info.base, /column,/align_center) ;plot base
  p_main = widget_base(info.base,/column,/align_center) ;parameter base
  obs_table = widget_base(p_main,/column,/align_center)
  lum_table = widget_base(p_main,/column,/align_center)
  dbase = widget_base(p_main,/column,/align_center) ; base for file dialogs
  button_base = widget_base(p_main,/row,/align_center) ; base for buttons

;Buttons at bottom 
  run_btn = widget_button(button_base,uvalue='go',value='Run Simulation',xsize=100,ysize=25)
  replot_btn = widget_button(button_base,uvalue='replot',value='Plot Last Run',xsize=100,ysize=25)
  quit_btn = widget_button(button_base,uvalue='quit',value='Quit',xsize=50,ysize=25)
  info_btn = widget_button(button_base,uvalue='info',value='Info',xsize=50,ysize=25)

;output file select
  obsname = fsc_fileselect(dbase,LabelName='Observation Save File',filename=ofile)
  oname = fsc_fileselect(dbase,LabelName='Output FITS file: ',filename=info.oname)

;survey parameter initialization
 if(file_test('params.save')) then begin
     restore,'params.save'
  endif else begin
     band1 = {wave:250.d0,fmin:25.0,ferr:6.2}
     band2 = {wave:350.d0,fmin:20.0,ferr:5.8}
     band3 = {wave:500.d0,fmin:15.0,ferr:6.2}
     bands = [band1,band2,band3]
  endelse

;luminosity function parameter initialization
  ldat={phi0:-2.2,lo:23.64,alpha:0.47,beta:2.88,p:-6.7,q:3.5}
  ldata=replicate(ldat,2) 
  ;the fixed values: =1 if held fixed, =0 if variable
  ldata(1).phi0=1
  ldata(1).lo=1
  ldata(1).alpha=1
  ldata(1).beta=1
  ldata(1).p=0
  ldata(1).q=0
  lrows=["Initial","Fixed"]

  bname = ["Band 1","Band 2","Band 3"]
  ocols = ["Wavelength (um)","Flux limit (mJy)","Standard Error (mJy)"]
  f = ['(e9.2)','(f7.4)','(f7.4)']
  fmt = [[f],[f],[f]]

;The Survey properties table
  lo = widget_label(obs_table,value="Survey Properties")
  ot = widget_table(obs_table,value=bands,column_labels=ocols,row_labels=bname,uvalue='ot',/editable,alignment=1,column_widths=[100,100,150,140],format=fmt)

;Luminosity Function Paramters
  l1 = widget_label(lum_table,value="Luminosity Function Parameters")
  t1 = widget_table(lum_table,value=ldata,column_labels=tag_names(ldat),row_labels=lrows,uvalue='t1',/editable,alignment=1,event_pro='bcheck')

;===================================================================================
;Show SED templates
;-----------------------------------------------------------------------------------
info.draw3 = widget_draw(info.base2, xsize=info.magnification*info.xsize2,$
                         ysize=info.magnification*info.ysize1 $
                         ,uvalue="DRAW_WINDOW3",retain=2 $
                         ,/button_events, keyboard_events=1,/tracking_events)

;initialize widget and establish control
widget_control,/realize,info.base,xoffset=0,yoffset=0
widget_control, info.draw3, get_value=win_id3 & info.win_id3=win_id3
xmanager,'simulation',info.base,/NO_BLOCK
wset,info.win_id3

;plot SEDs (generalize to any SED template file
templ=mrdfits('sf_templates.fits')
plot,templ[*,0],templ[*,1],/xlog,/ylog,yrange=[1.d20,1.d28],ystyle=1,xtitle=TeXtoIDL('\lambda [\mum]'),ytitle=TeXtoIDL('L_{\nu} [W/Hz]')
for ipl=1,13 do oplot,templ[*,0],templ[*,ipl+1]

END

;widget event handling routine
PRO simulation_event,ev
  COMMON simulation_com,cdir,info,ot,ldist,lum,settings,bands

  ; get event identifier
  widget_control,ev.id,get_uvalue=uvalue

  CASE uvalue OF
     'go'  : begin ;save settings, intialize FITS file, and pass to C++ fitting routine
        save,ldist,lum,settings,bands,filename='params.save' ;save parameters
        
        ;widget_control,ot,get_value=value ;get observation table values, print, value
        ;widget_control,t1,get_value=value ;get lumfunct table values

        ;make observation FITS file
        widget_control,ot,get_value=bvals
        wave = [bvals[0].wave,bvals[1].wave,bvals[2].wave]
        flux_min = [bvals[0].fmin,bvals[1].fmin,bvals[2].fmin]
        flux_err = [bvals[0].ferr,bvals[1].ferr,bvals[2].ferr]

        widget_control,obsname,get_value=value
        print,value
        if (file_test(info.ofile)) then begin
           restore,info.ofile
           values = dblarr(n_elements(wave),n_elements(f1))
        endif else begin
           print,'Error: file "observation.save" not found'
           print,'Simulation attempt unsuccessful'
           break
        endelse

        for i=0,n_elements(f1)-1 do begin
           values(0,i) = f1[i]
           values(1,i) = f2[i]
           values(2,i) = f3[i]
        endfor
        
        sxaddpar, hdr, 'DATE', systime(),'Date of creation'
        sxaddpar, hdr, 'FLUX_NUM',n_elements(f1),'Number of observations included in file'
        sxaddpar, hdr, 'WAVE_1',wave[0],'Wavelength corresponding to first flux'
        sxaddpar, hdr, 'WAVE_2',wave[1],'Wavelength corresponding to second flux'
        sxaddpar, hdr, 'WAVE_3',wave[2],'Wavelength corresponding to third flux'
        sxaddpar, hdr, 'W1_FMIN',flux_min[0],'Wavelength corresponding to first flux'
        sxaddpar, hdr, 'W2_FMIN',flux_min[1],'Wavelength corresponding to second flux'
        sxaddpar, hdr, 'W3_FMIN',flux_min[2],'Wavelength corresponding to third flux'
        sxaddpar, hdr, 'W1_FERR',flux_err[0],'Wavelength corresponding to first flux'
        sxaddpar, hdr, 'W2_FERR',flux_err[1],'Wavelength corresponding to second flux'
        sxaddpar, hdr, 'W3_FERR',flux_err[2],'Wavelength corresponding to third flux'
                
        mwrfits,values,'observation.fits',hdr,/create

        ;make model fits file

        widget_control,t1a,get_value=zvals
        widget_control,t1b,get_value=pvals

        restore,info.mfile
        tmp_info = size(templates)
        mod_nums = n_elements(templates)/n_elements(wave)
        templates = reform(templates,n_elements(templates))
        
        sxaddpar, hdr2, 'DATE', systime(),'Date of creation'
        sxaddpar, hdr2, 'P_NUM',info.pnum,'Number of Model Parameters'
        sxaddpar, hdr2, 'Z_NUM',info.znum,'Number of Redshift Distributions'

        n_els = n_elements(wave)
        wave_min = wave[0]
        wave_max = wave[n_els-1]
        wave_sep = (wave_max-wave_min)/(n_els-1)
        
        sxaddpar, hdr2, 'WAVE_MIN',wave_min,'Domain Lower Bound'
        sxaddpar, hdr2, 'WAVE_MAX',wave_max,'Domain Upper Bound'
        sxaddpar, hdr2, 'WAVE_SEP',wave_sep,'Domain Step Size'
        
        for i=0,info.znum-1 do begin
           z_tmp = strcompress(string(i),/remove_all)
           sxaddpar, hdr2, 'Z'+z_tmp+'_MAX',zvals[i].max,'Parameter '+z_tmp+' Upper Limit'
           sxaddpar, hdr2, 'Z'+z_tmp+'_MIN',zvals[i].min,'Parameter '+z_tmp+' Lower Limit'
           sxaddpar, hdr2, 'Z'+z_tmp+'_MEAN',zvals[i].mean,'Parameter '+z_tmp+' Distribution Mean'
           sxaddpar, hdr2, 'Z'+z_tmp+'_SIGMA',zvals[i].sigma,'Parameter '+z_tmp+' Distribution Width'
           sxaddpar, hdr2, 'Z'+z_tmp+'_FIXED',zvals[i].fixed,'Parameter '+z_tmp+' Fixed or Variable'
        endfor

        for i=0,info.pnum-1 do begin
           p_tmp = strcompress(string(i),/remove_all)
           sxaddpar, hdr2, 'P'+p_tmp+'_SIZE',tmp_info[i+1],'Parameter '+p_tmp+' Upper Limit'
           sxaddpar, hdr2, 'P'+p_tmp+'_MAX',pvals[i].max,'Parameter '+p_tmp+' Upper Limit'
           sxaddpar, hdr2, 'P'+p_tmp+'_MIN',pvals[i].min,'Parameter '+p_tmp+' Lower Limit'
           sxaddpar, hdr2, 'P'+p_tmp+'_MEAN',pvals[i].mean,'Parameter '+p_tmp+' Distribution Mean'
           sxaddpar, hdr2, 'P'+p_tmp+'_SIGMA',pvals[i].sigma,'Parameter '+p_tmp+' Distribution Width'
           sxaddpar, hdr2, 'P'+p_tmp+'_FIXED',pvals[i].fixed,'Parameter '+p_tmp+' Fixed or Variable'
        endfor
        
        widget_control,t1,get_value=lparam

        sxaddpar,hdr2,'PHI0',lparam[0],'Luminosity Function Normalization'
        sxaddpar,hdr2,'L0',lparam[1],'Luminosity Function Knee'
        sxaddpar,hdr2,'ALPHA',lparam[2],'Luminosity Function upper slope'
        sxaddpar,hdr2,'BETA',lparam[3],'Luminosity Function lower slope'
        sxaddpar,hdr2,'P',lparam[4],'Luminosity Function PHI evolution term'
        sxaddpar,hdr2,'Q',lparam[5],'Luminosity Function L evolution term'

        print,"Luminosity function parameters: ",lparam

        mwrfits,templates,'model.fits',hdr2,/create

;Actual simulation
        
        spawn,cdir+'z_fit'

        read_output
        graphs
     end
     'replot': begin
        read_output
        graphs
     end
     'quit': begin
        widget_control,ev.top,/destroy
     end
     't1' :
     'ot': widget_control,ot,get_value=bands
     'size': widget_control,size,get_value=settings
     'info': fit_info
     ELSE:
  ENDCASE
END

pro bcheck,ev

  i = ev.y

  widget_control,ev.id,get_value=value
  if(value[i].fixed lt 0) then begin
     value[i].fixed=0
     widget_control,ev.id,set_value=value
  endif else if(value[i].fixed gt 1) then begin
     value[i].fixed=1
     widget_control,ev.id,set_value=value
  endif

end

pro fit_info
  
  main = widget_base(title="Information about Fitting Methodology",/column)
  t1 = widget_text(main,value="this is text")
  bt = widget_button(main,uvalue='close',value='Close',xsize=50,ysize=25)

  widget_control,main,/realize
  xmanager,'fit_info',main,/no_block

end

pro fit_info_event,ev

  widget_control,ev.id,get_uvalue=uvalue

  case uvalue of
     'close': begin
        widget_control,ev.top,/destroy
     end
  endcase

end

pro read_output
  
  !p.thick=5
  !x.thick=5
  !y.thick=5
  !p.charthick=5
  !p.charsize=1.5

  comp = mrdfits('output.fits',0,head,/silent)
  model = mrdfits('output.fits',1,/silent)
  obs = mrdfits('output.fits',2,/silent)

  xysize = fxpar(head,'DIM')
  hist_min = fxpar(head,'H_MIN')
  hist_max = fxpar(head,'H_MAX')
  binsize = fxpar(head,'BINSIZE')

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
  resmod[neg] = 0
  resobs[pos] = 0

  plot,[hist_min,hist_max],[hist_min,hist_max],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}')
  hb = 0.5

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

end

pro graphs

  !p.thick=0
  !x.thick=0
  !y.thick=0
  !p.charthick=0
  !p.charsize=0

  size_screen=get_screen_size()
  size_screen=size_screen*0.8
  gmain = widget_base(title='Simulation Output',/column,xsize=size_screen[0],ysize=size_screen[1])
  r1 = widget_base(gmain,/row) ;,xsize=size_screen[0],ysize=size_screen[1])
  r2 = widget_base(gmain,/row) ;,xsize=size_screen[0],ysize=size_screen[1])
  r2b = widget_base(gmain,/row) ;,xsize=size_screen[0],ysize=size_screen[1])
  r3 = widget_base(gmain,/row) ;,xsize=size_screen[0],ysize=size_screen[1])

  widget_control,gmain,set_uvalue=mnum
  xdim = fix(size_screen[0]/3.0)
  ydim = fix(size_screen[1]/3.0)

  lumfunct = widget_draw(r1,xsize=xdim,ysize=ydim)
  redshift = widget_draw(r1,xsize=xdim,ysize=ydim)
  models = widget_draw(r1,xsize=xdim,ysize=ydim)
  dcount1 = widget_draw(r2,xsize=xdim,ysize=ydim)
  dcount2 = widget_draw(r2,xsize=xdim,ysize=ydim)
  dcount3 = widget_draw(r2,xsize=xdim,ysize=ydim)
  sim_colors = widget_draw(r2b,xsize=xdim,ysize=ydim)
  obs_colors = widget_draw(r2b,xsize=xdim,ysize=ydim)
  comp_colors = widget_draw(r2b,xsize=xdim,ysize=ydim)

  refresh = widget_button(r3,uvalue='refresh',value='Refresh')
  close = widget_button(r3,uvalue='close',value='Close')
  quit = widget_button(r3,uvalue='quit',value='Quit')

  widget_control,gmain,/realize
  xmanager,'graphs',gmain,/no_block

  set_plot,'x'
  loadct,0,/silent

  dists = mrdfits('output.fits',3,head,/silent)

  gpts = where(dists.f3 gt 0)
  f1 = dists[gpts].f1
  f2 = dists[gpts].f2
  f3 = dists[gpts].f3
  z = dists[gpts].z
  m = dists[gpts].m
  lum = dists[gpts].lum

  pnum_out = fxpar(head,'tfields')

  for i=6,pnum_out-1 do begin
     if (i eq 6) then begin
        hists = histogram(dists[gpts].(i),binsize=0.1,locations=xh,min=0,max=10)
        hmax = max(hists)
        xmax = max(xh)
     endif else begin
        h = histogram(dists[gpts].(i),binsize=0.1,locations=xh,min=0)
        if(max(xh) gt xmax) then begin
           xmax = max(xh)
           for j=6,i-1 do begin
              if(i eq 6) then begin
                 hists = histogram(dists[gpts].(i),binsize=0.1,locations=xh,min=0,max=xmax)
              endif else begin
                 htemp = histogram(dists[gpts].(i),binsize=0.1,locations=xh,min=0,max=xmax)
                 hists = [hists,htemp]
              endelse
           endfor
        endif else begin
           h = histogram(dists[gpts].(i),binsize=0.1,locations=xh,min=0,max=xmax)
        endelse
        temp = max(h)
        hmax = (temp gt hmax) ? temp : hmax
        hists = [[hists],[h]]
     endelse
  endfor     

  widget_control,lumfunct,get_value=index
  wset,index
  h = histogram(alog10(lum),nbins=50,locations=xh)
  plot,xh,h,psym=10,xstyle=1,yrange=[0.1,1000],ystyle=1,xtitle='Log(Luminosity (W/Hz))',ytitle='dN/(dL/dHz)',/ylog,title='Luminosity Function'

  widget_control,redshift,get_value=index
  wset,index
  h = histogram(z,binsize=0.1,locations=xh,min=0,max=10)
  plot,xh,h,psym=10,xrange=[0,10],xstyle=1,xtitle='z',ytitle='dN/dz',title='Redshift Distribution'

  widget_control,models,get_value=index
  wset,index

  help,hists
  for i=6,pnum_out-1 do begin
     if(i eq 6) then begin
        plot,xh,hists(*,i-6),psym=10,xstyle=1,xtitle='m',ytitle='dN/dm',xrange=[0,xmax],title='Model Distribution',yrange=[0,hmax*1.2]
     endif else begin
        oplot,xh,hists(*,i-6),linestyle=(i-6),psym=10
     endelse
  endfor

  widget_control,dcount1,get_value=index
  wset,index
;Herschel ATLAS counts at 250,350 and 500 (Clements et al. 2010)
  readcol,'counts_clements10.dat',skipline=2,numline=16,flux,nbin,corr,int_counts,int_err,diff_counts,diff_err

  flux=flux/1.d3
  plot,flux,diff_counts,psym=sym(1),symsize=2,xtitle=TeXtoIDL('F_{250}[Jy]'),ytitle=TeXtoIDL('(dN/dS)S^{2.5} [gal ster^{-1} J^{1.5}]'),/xlog,/ylog,title='Band 1 Counts',yrange=[5.d2,1.d5],ystyle=1

  oploterr,flux,diff_counts,diff_err
  h = histogram(alog10(f1),nbins=50,locations=xh,min=-3,max=1)
  ;I think the area covered here is not
;quite right needs to be sorted out better, but for now just scale up
  h=h*2.d3
  ;dlogS=dS/S
  df1=(shift(xh,-1)-xh)*10.^(xh)
  dcounts=(h/df1)*f1^(2.5)
  xh=10.^(xh)
  oplot,xh,dcounts,psym=10

  widget_control,dcount2,get_value=index
  wset,index

  h = histogram(alog10(f2),nbins=50,locations=xh,min=-3,max=0)
  pts = where(h le 0)
  h[pts] = 0.01
  plot,xh,h,psym=10,/ylog,xrange=[-3,0],yrange=[1e-1,1e3],ystyle=1,xstyle=1,xtitle='Flux [Log(Jy)]',ytitle='dN/dS (Log)',title='Band 2 Counts'

  widget_control,dcount3,get_value=index
  wset,index

  h = histogram(alog10(f3),nbins=50,locations=xh,min=-3,max=0)
  pts = where(h le 0)
  h[pts] = 0.01
  plot,xh,h,psym=10,/ylog,xrange=[-3,0],yrange=[1e-1,1e3],ystyle=1,xstyle=1,xtitle='Flux [Log(Jy)]',ytitle='dN/dS (Log)',title='Band 3 Counts'

  comp = mrdfits(info.oname,0,head,/silent)
  model = mrdfits(info.oname,1,/silent)
  obs = mrdfits(info.oname,2,/silent)

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

        if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
     endfor
     oplot,[a[i+1],a[i+1]],[hist_min,hist_max],linestyle=1
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

        if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
     endfor
     oplot,[a[i+1],a[i+1]],[hist_min,hist_max],linestyle=1
  endfor

  widget_control,comp_colors,get_value=index
  wset,index
    
  pos = where(comp gt 0)
  neg = where(comp lt 0)
  resmod = comp
  resobs = abs(comp)
  resmod[neg] = 0
  resobs[pos] = 0

  plot,[hist_min,hist_max],[hist_min,hist_max],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}'),title='Color Distribution Comparison'
  hb = 0.5

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

end

pro graphs_event,ev
  COMMON simulation_com,cdir,info,ot,ldist,lum,settings,bands

  widget_control,ev.id,get_uvalue=uvalue

  case uvalue of
     'refresh': begin
        widget_control,ev.top,get_uvalue=mnum
        widget_control,ev.top,/destroy
        graphs
     end
     'close': widget_control,ev.top,/destroy
     'quit': begin
        widget_control,ev.top,/destroy
        widget_control,info.base,/destroy
     end
  endcase

end
