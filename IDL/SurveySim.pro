;==================================================================
; Writen by Noah Kurinsky, version recent as of 8/21/14
;==================================================================

pro SurveySim
  COMMON simulation_com,info,parameters
   
  CD, Current=thisdir 
  plot_settings                 ;will reset the global settings
  
;==================================================================
;the INFO structure holds the key widget control parameters as well as 
;the basic simulation settings
;------------------------------------------------------------------
  info = make_info_struct()
  
;====================================================================
; Colors
;--------------------------------------------------------------------

  loadct,3,/silent
  tvlct,red,green,blue,/get
  info.ncolors=!d.TABLE_SIZE
  red[info.ncolors-1]=0B
  green[info.ncolors-1]=255B
  blue[info.ncolors-1]=0B
  tvlct,red,green,blue
  
; Screen size
  size_screen = get_screen_size()
  size_screen=size_screen*0.8
  spectrum_xsize=size_screen(0)*2./3. & spectrum_ysize=size_screen(1)/5.
  
  plotSize = 100.
  info.magnification = size_screen[1]/(plotSize+info.ysize1)*0.9
  
  ;widget base initialization
  ;main base
  info.base = widget_base(title='Model Setup and Initialization',mbar=mbar,/column,/align_center) 
  ;plot base
  info.base2= widget_base(info.base, /column,/align_center)
  ;parameter base
  info.p_main = widget_base(info.base,/column,/align_center)
  
  info.obs_table = widget_base(info.p_main,/column,/align_center)
  info.lum_table = widget_base(info.p_main,/column,/align_center)
  info.sim_table = widget_base(info.p_main,/column,/align_center)
  info.dbase = widget_base(info.p_main,/column,/align_center) ; base for file dialogs
  info.button_base = widget_button(mbar,/menu,value="Menu")   ; base for buttons
  
  ;=============================================================
  ;survey parameter initialization
  ;-------------------------------------------------------------
  ;attempt to get saved parameters, create new if not saved
  parameters = load_parameters(thisdir+"/params.save")
  filter_names = filter_list("/usr/local/surveysim/filters/filterlib.txt")
  filter_names = ["Select Filter",filter_names]
  
  ;The filter properties table
  lo = widget_label(info.obs_table,value="Survey Properties")
  obs_row = widget_base(info.obs_table,/row)
  filter_dialogs = widget_base(obs_row,/column,/align_center)
  info.fd1 = widget_combobox(filter_dialogs, $
                             value=filter_names,$
                             title="Filter 1",$
                             uvalue="fd1")
  info.fd2 = widget_combobox(filter_dialogs, $
                             value=filter_names,$
                             title="Filter 2",$
                             uvalue="fd2")
  info.fd3 = widget_combobox(filter_dialogs, $
                             value=filter_names,$
                             title="Filter 3",$
                             uvalue="fd3")
  info.ot = widget_table(obs_row,$
                         uvalue='ot',/editable,alignment=1,$
                         value=parameters.filters.properties,$
                         column_labels=["Flux limit (mJy)","Standard Error (mJy)"],$
                         row_labels=["Filter 1","Filter 2","Filter 3"],$
                         column_widths=[100,100],$
                         format=[[f],[f]],$
                         scr_xsize=425,scr_ysize=95)

  widget_control, info.fd1, set_combobox_select=parameters.filters[0].filter_id
  widget_control, info.fd2, set_combobox_select=parameters.filters[1].filter_id
  widget_control, info.fd3, set_combobox_select=parameters.filters[2].filter_id

  ;Survey Model Parameters
  l1 = widget_label(info.lum_table,value="Survey Model Parameters")
  info.t1 = widget_table(info.lum_table, $
                         value=parameters.lumpars.pars,$
                         row_labels=parameters.lumpars.name,$
                         row_labels=tag_names(parameters.lumpars.pars),$
                         uvalue='t1',$
                         /editable,alignment=1,$
                         format=make_array(shape(parameters.lumpars.pars),'(f5.2)'),$
                         scr_xsize=537,scr_ysize=132)
  
  ;Simulation Parameters
  l2 = widget_label(info.sim_table,value="Simulation Settings")
  info.t2 = widget_table(info.sim_table,$ 
                         value=parameters.surveyData,$
                         column_labels=["Area (sdeg)","Z Min","Z Max","Z Binsize","Run Number"],$
                         /no_row_headers,uvalue='t2',/editable,alignment=1,$
                         format=['(f5.2)','(f5.2)','(f5.2)','(f5.2)','(e9.2)'],$
                         column_widths=[100,100,100,100,100],$
                         scr_xsize=505,scr_ysize=55)
  
  
;===========================================================================
;Buttons in Menu
;---------------------------------------------------------------------------
  sim_opt= widget_button(mbar,/menu,value="Simulation")
  run_btn = widget_button(sim_opt,uvalue='go',value='Run')
  set_btn = widget_button(sim_opt,uvalue='settings',value='Settings')
  plots = widget_button(mbar,/menu,value="Plots")
  replot_btn = widget_button(plots,uvalue='replot',value='Simulation Output')
  diag_btn = widget_button(plots,uvalue='diag',value='MCMC Diagnostics')
  info_btn = widget_button(info.button_base,uvalue='info',value='About')
  save_btn = widget_button(info.button_base,uvalue='save',value='Save Settings')
  quit_btn = widget_button(info.button_base,uvalue='quit',value='Quit')
  
;============================================================================
;Show SED templates
;----------------------------------------------------------------------------
  info.draw3 = widget_draw(info.base2, xsize=info.magnification*info.xsize2,$
                           ysize=info.magnification*info.ysize1, $
                           uvalue="DRAW_WINDOW3",retain=2, $
                           /button_events, keyboard_events=1,/tracking_events)
  
  ;initialize widget and establish control
  widget_control,/realize,info.base,xoffset=0,yoffset=0,Set_UValue=theObject
  widget_control, info.draw3, get_value= info.win_id3 
  xmanager,'SurveySim',info.base,/NO_BLOCK
  wset,info.win_id3

  ;get new sedfile if file does not exist (option to ignore too)
  if not file_test(files.sedfile) then begin
     file.sedfile = ask_for_file("SED File "+files.sedfile+" not found, please enter new SED file",files.sedfile)
  endif
  
  ;plot SEDs
  set_plot,'x'
  device,decomposed=0
  plot_settings,plot_type='x'
  if file_test(files.sedfile) then begin

     templ=mrdfits(files.sedfile,/silent) 
     loadct,1,/silent

     plot,templ[*,0],templ[*,1],/xlog,/ylog,$
          yrange=[1.d20,1.d28],ystyle=1,ytitle=TeXtoIDL('L_{\nu} [W/Hz]'),$
          xtitle=TeXtoIDL('\lambda [\mum]')

     for ipl=1,13 do oplot,templ[*,0],templ[*,ipl+1]
  endif else begin
     print,'File '+files.sedfile+' does not exist, skipping plot'
  endelse
  
END

;widget event handling routine
PRO SurveySim_event,ev
  COMMON simulation_com

  ; get event identifier
  widget_control,ev.id,get_uvalue=uvalue

  CASE uvalue OF
     'save' : save,ldata,bands,sdat,cdat,msettings,files,filename='params.save' ;save parameters
     'go'  : begin              ;save settings, intialize FITS file, and pass to C++ fitting routine
        save,ldata,bands,sdat,cdat,msettings,files,filename='params.save' ;save parameters
        
                                ;update observation FITS file
        widget_control,info.ot,get_value=bvals
        flux_min = [bvals[0].fmin,bvals[1].fmin,bvals[2].fmin]
        flux_err = [bvals[0].ferr,bvals[1].ferr,bvals[2].ferr]

        hdr = headfits(files.ofile)
        hnum = sxpar(hdr,"FHDU")
        hdr = headfits(files.ofile,exten=hnum)

        ;add filters
        sxaddpar, hdr, 'F1MIN',flux_min[0],'Flux cutoff, first column'
        sxaddpar, hdr, 'F2MIN',flux_min[1],'Flux cutoff, second column'
        sxaddpar, hdr, 'F3MIN',flux_min[2],'Flux cutoff, third column'
      
        modfits,files.ofile,0,hdr,exten_no=hnum

        ;make model fits file

        widget_control,info.t1,get_value=lparam
        pars = lparam(0)
        fixed = lparam(1)
        min = lparam(2)
        max = lparam(3)
        widget_control,info.t2,get_value=sparam
        sdat = sparam
        widget_control,info.t3,get_value=cparam
        cdat = cparam

        sxaddpar,hdr2,'DATE',systime(),'Date of creation'

        sxaddpar,hdr2,'PHI0',pars.phi0,'Luminosity Function Normalization'
        sxaddpar,hdr2,'PHI0_FIX',fixed.phi0,'Fix Phi0 (Y=1/N=0)'
        sxaddpar,hdr2,'PHI0_MIN',min.phi0,'Minimum Phi0 value'
        sxaddpar,hdr2,'PHI0_MAX',max.phi0,'Maximum Phi0 value'

        sxaddpar,hdr2,'L0',pars.lo,'Luminosity Function Knee'
        sxaddpar,hdr2,'L0_FIX',fixed.lo,'Fix L0 (Y=1/N=0)'
        sxaddpar,hdr2,'L0_MIN',min.lo,'Minimum L0 value'
        sxaddpar,hdr2,'L0_MAX',max.lo,'Maximum L0 value'

        sxaddpar,hdr2,'ALPHA',pars.alpha,'Luminosity Function upper slope'
        sxaddpar,hdr2,'ALPHA_FIX',fixed.alpha,'Fix Alpha (Y=1/N=0)'
        sxaddpar,hdr2,'ALPHA_MIN',min.alpha,'Minimum Alpha value'
        sxaddpar,hdr2,'ALPHA_MAX',max.alpha,'Maximum Alpha value'

        sxaddpar,hdr2,'BETA',pars.beta,'Luminosity Function lower slope'
        sxaddpar,hdr2,'BETA_FIX',fixed.beta,'Fix Beta (Y=1/N=0)'
        sxaddpar,hdr2,'BETA_MIN',min.beta,'Minimum Beta value'
        sxaddpar,hdr2,'BETA_MAX',max.beta,'Maximum Beta value'

        sxaddpar,hdr2,'P',pars.p,'Luminosity Function PHI evolution term'
        sxaddpar,hdr2,'P_FIX',fixed.p,'Fix P (Y=1/N=0)'
        sxaddpar,hdr2,'P_MIN',min.p,'Minimum P value'
        sxaddpar,hdr2,'P_MAX',max.p,'Maximum P value'

        sxaddpar,hdr2,'Q',pars.q,'Luminosity Function L evolution term'
        sxaddpar,hdr2,'Q_FIX',fixed.q,'Fix Q (Y=1/N=0)'
        sxaddpar,hdr2,'Q_MIN',min.q,'Minimum Q value'
        sxaddpar,hdr2,'Q_MAX',max.q,'Maximum Q value'

        sxaddpar,hdr2,'ZCUT',pars.zcut,'Luminosity Function z evolution limit'
        sxaddpar,hdr2,'ZCUT_FIX',fixed.zcut,'Fix ZCUT (Y=1/N=0)'
        sxaddpar,hdr2,'ZCUT_MIN',min.zcut,'Minimum ZCUT value'
        sxaddpar,hdr2,'ZCUT_MAX',max.zcut,'Maximum ZCUT value'

        sxaddpar,hdr2,'CEXP',cdat.a0,'Intrinsic luminosity evolution term'
        sxaddpar,hdr2,'CEXP_FIX',cdat.fixed,'Fix CEXP (Y=1/N=0)'
        sxaddpar,hdr2,'CEXP_MIN',cdat.amin,'Minimum CEXP value'
        sxaddpar,hdr2,'CEXP_MAX',cdat.amax,'Maximum CEXP value'

        sxaddpar,hdr2,'RUNS',sdat.runs,'Number of Runs'
        sxaddpar,hdr2,'ZMIN',sdat.zmin,'Minimum Redshift Value'
        sxaddpar,hdr2,'ZMAX',sdat.zmax,'Maximum Redshift Value'
        sxaddpar,hdr2,'DZ',sdat.dz,'Redshit Bin Width'
        sxaddpar,hdr2,'AREA',sdat.area,'Observed Solid Angle'

        sxaddpar,hdr2,'NCHAIN',msettings.nchain,'Chain Number'
        sxaddpar,hdr2,'TMAX',msettings.tmax,'Starting Anneal Temperature'
        sxaddpar,hdr2,'ANN_PCT',msettings.acceptpct,'Ideal Acceptance Percentage'
        sxaddpar,hdr2,'CONV_CONF',msettings.conv_conf,'Convergence CI Setting'
        sxaddpar,hdr2,'CONV_RMAX',msettings.conv_rmax,'Convergence Rmax Criterion'
        sxaddpar,hdr2,'CONV_STEP',msettings.conv_step,'Iterations between convergence checks'
        sxaddpar,hdr2,'BURN_STEP',msettings.burn_step,'Iterations between anneal calls in burn-in'
        sxaddpar,hdr2,'BURNVRUN',msettings.burn_ratio,'Ratio of normal to burn-in steps'
        sxaddpar,hdr2,'ANN_RNG',msettings.pct_range,'Range within which to maintain acceptance, from ideal'
        sxaddpar,hdr2,'PRINT',info.print,'Whether to Print Debug MSGs'

        templates = [0]
        mwrfits,templates,files.mfile,hdr2,/create

;===============================================================
;Run the actual simulation
;---------------------------------------------------------------
        args = files.ofile+' '+files.mfile+' '+files.sedfile+' '+files.oname
        spawn,'fitter '+args

        read_output
        graphs
     end
     'settings': begin
        settings
        
        wset,info.win_id3
        if file_test(files.sedfile) then begin
           templ=mrdfits(files.sedfile,/silent) 
           loadct,1,/silent
           plot,templ[*,0],templ[*,1],/xlog,/ylog,yrange=[1.d20,1.d28],ystyle=1,xtitle=TeXtoIDL('\lambda [\mum]'),ytitle=TeXtoIDL('L_{\nu} [W/Hz]')
           for ipl=1,13 do oplot,templ[*,0],templ[*,ipl+1]
        endif else begin
           print,'File '+files.sedfile+' does not exist, skipping plot'
        endelse

     end
     'diag': diagnostics
     'replot': begin
        read_output
        graphs
     end
     'quit': begin
        widget_control,ev.top,/destroy
     end
     'ot': widget_control,info.ot,get_value=bands
     't1': begin
        widget_control,info.t1,get_value=ldata
        
        i = ev.y
        j = ev.x        

        if (i eq 1) then begin
           if(ldata[i].(j) lt 1) then begin
              ldata[i].(j)=0
              widget_control,ev.id,set_value=ldata
           endif else if(ldata[i].(j) gt 1) then begin
              ldata[i].(j)=1
              widget_control,ev.id,set_value=ldata
           endif
        endif

     end
     't2': widget_control,info.t2,get_value=sdat
     't3': widget_control,info.t3,get_value=cdat
     'info': SurveySim_info
     ELSE:
  ENDCASE

END

pro read_output
  
  COMMON simulation_com

  print,files.oname
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

  count_dists = mrdfits(files.oname,6,head,/silent)

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

        ;if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
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

        ;if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
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
           
           ;if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
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
  
  device,/close

  if(file_test('plots')) then begin
     file_delete,'plots',/recursive
  endif
  file_mkdir,'plots'
  file_move,'*.eps','plots'

end

pro graphs

  COMMON simulation_com

  size_screen=get_screen_size()
  size_screen_alt = size_screen*0.85
  size_screen = size_screen*0.8
  gmain = widget_base(title='Simulation Output',/column,xsize=size_screen[0],ysize=size_screen_alt[1])
  r1 = widget_base(gmain,/row) ;,xsize=size_screen[0],ysize=size_screen[1])
  r2 = widget_base(gmain,/row) ;,xsize=size_screen[0],ysize=size_screen[1])
  r2b = widget_base(gmain,/row) ;,xsize=size_screen[0],ysize=size_screen[1])
  r3 = widget_base(gmain,/row) ;,xsize=size_screen[0],ysize=size_screen[1])

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
  xmanager,'graphs',gmain,/no_block
  
  set_plot,'x'
  device,decomposed=0
  plot_settings,plot_type='x'
  loadct,0,/silent

  res = mrdfits(files.oname,0,head,/silent)

  alpha = fxpar(head,'ALPHA')
  beta = fxpar(head,'BETA')
  phi0 = fxpar(head,'PHI0')
  L0 = fxpar(head,'L0')
  p = fxpar(head,'P')
  q = fxpar(head,'Q')

  omega_m = 0.28
  omega_l = 0.72
  
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

  count_dists = mrdfits(files.oname,6,head,/silent)

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

  pnum_out = fxpar(head,'tfields')

  widget_control,lumfunct,get_value=index
  wset,index
  device,decomposed=0

  plot,[8,13],[1e-10,1e0],/ylog,/nodata,xtitle=textoidl('log_{10}(L_{fir}) [L_{sun}]'),ytitle=textoidl('log_{10}(N(L_{fir})/\Omega dV_C)'),ystyle=1
  lums = indgen(21)/4.0+8.0

  dz = 0.5

  loadct,39,/silent
  for zi=dz,5.0,dz do begin
     if zi le 2.0 then begin
        t1 = (10^phi0)*((1+zi)^p)
        t2 = (10^L0)*((1+zi)^q)
     endif else begin
        t1 = (10^phi0)*((3.0)^p)
        t2 = (10^L0)*((3.0)^q)
     endelse

     r = 10^lums/t2
     nsrcs=t1/(r^alpha+r^beta)

     oplot,lums,nsrcs,color=zi*40+20
  endfor

  for i=0.0,5.0,0.1 do begin
     ind = (9+2*i/5)
     oplot,[ind,ind+0.02,ind,ind+0.02,ind,ind+0.02],[10^(-8),10^(-8),10^(-7.9),10^(-7.9),10^(-7.8),10^(-7.8)],psym=2,color=i*40+20,symsize=0.75
  endfor
  
  loadct,0,/silent
  oplot,[8,13],[1,1],linestyle=1

  xyouts,9.0,10^(-7.2),'Redshift'
  xyouts,9.0,10^(-7.6),'0'
  xyouts,10.9,10^(-7.6),'5'

  widget_control,redshift,get_value=index
  wset,index
  h = histogram(z,binsize=0.1,locations=xh,min=0.2,max=5.0)
  plot,xh,h,psym=10,xrange=[0,max(z)],xstyle=1,xtitle='z',ytitle='dN/dz',title='Redshift Distribution'

  widget_control,models,get_value=index
  wset,index

  f = alog10(lum)
  h = histogram(f,nbins=20,locations=xh,min=min(f),max=max(f))
  binsize=xh[1]-xh[0]
  xpts=10.0^(xh+binsize/2)
  plot,xpts,h,psym=2,/xlog,/ylog,xstyle=1,ystyle=0,title="Luminosity Distribution",xtitle="L",ytitle="dN/dL"

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

        ;if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
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

        ;if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
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
           
           ;if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
        endfor
        ;oplot,[a[i+1],a[i+1]],[hist_min,hist_max],linestyle=1
     endfor
  endfor

end

pro graphs_event,ev
  COMMON simulation_com

  widget_control,ev.id,get_uvalue=uvalue

  case uvalue of
     'diags' : begin
        widget_control,ev.top,/destroy
        diagnostics
     end
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


pro diagnostics

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
  xmanager,'diagnostics',dmain,/no_block
  
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

pro diagnostics_event,ev
  COMMON simulation_com

  widget_control,ev.id,get_uvalue=uvalue

  case uvalue of
     'graphs': begin
        widget_control,ev.top,/destroy
        graphs
     end
     'refresh': begin
        widget_control,ev.top,get_uvalue=mnum
        widget_control,ev.top,/destroy
        diagnostics
     end
     'close': widget_control,ev.top,/destroy
     'quit': begin
        widget_control,ev.top,/destroy
        widget_control,info.base,/destroy
     end
  endcase

end

