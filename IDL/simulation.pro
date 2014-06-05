
;==================================================================
; Writen by Noah Kurinsky, version recent as of 10/17/13
; Edited by Anna Sajina, December 2012
; certain files in the same directory are required for proper
; function of this widget:
;
;==================================================================

pro simulation

  COMMON simulation_com,info,ldata,ldat0,sdat,cdat,bands,msettings,files
  plot_settings ;will reset the global settings

;==================================================================
;the INFO structure holds the key widget control parameters as well as 
;the basic simulation settings
;------------------------------------------------------------------
  info={$                       ;widget settings
       base: 0L, $
       base1: 0L, $
       base2: 0L, $
       base_out: 0L, $
       base1_out: 0L, $
       base2_out: 0L, $
       xsize1:150., $
       ysize1:100., $
       xsize2:150., $
       ysize2:100., $
       magnification:1.00, $  
       draw1:0L, $
       draw2:0L, $
       draw3:0L, $
       win_id1:0L, $
       win_id2:0L, $
       win_id3:0L, $
       ncolors:0L, $
;widget bases
       p_main:0L, $
       obs_table:0L, $
       lum_table:0L, $
       sed_table:0L, $
       sim_table:0L, $
       dbase:0L, $
       button_base:0L, $
;input tables and files
       obsname:0L, $
       sfile:0L, $
       oname:0L, $
       mname:0L, $
       ot:0L, $
       t1:0L, $
       t2:0L, $
       t3:0L, $
       tset:0L, $
       dprint:0L, $
       print:1}
  
;====================================================================
; Colors
;--------------------------------------------------------------------
;device, truecolor, depth=24 ;decompose=0
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
  info.base = widget_base(title='Model Setup and Initialization',mbar=mbar,/column,/align_center) ;main base
  info.base2= widget_base(info.base, /column,/align_center)                                       ;plot base
  info.p_main = widget_base(info.base,/column,/align_center)                                      ;parameter base
  info.obs_table = widget_base(info.p_main,/column,/align_center)
  info.lum_table = widget_base(info.p_main,/column,/align_center)
  info.sed_table = widget_base(info.p_main,/column,/align_center)
  info.sim_table = widget_base(info.p_main,/column,/align_center)
  info.dbase = widget_base(info.p_main,/column,/align_center) ; base for file dialogs
  info.button_base = widget_button(mbar,/menu,value="Menu")   ; base for buttons
  
;=============================================================
;survey parameter initialization
;-------------------------------------------------------------
;don't want to use the existing params.save as the wavelength
;is in meters, whereas we use microns normally
  if(file_test('params.save')) then begin
     restore,'params.save'
  endif else begin
     band1 = {wave:250.d0,fmin:25.0,ferr:6.2}
     band2 = {wave:350.d0,fmin:20.0,ferr:5.8}
     band3 = {wave:500.d0,fmin:15.0,ferr:6.2}
     bands = [band1,band2,band3]
     
                                ;luminosity function parameter initialization
     ldat0={phi0:-2.2,lo:10.14,alpha:0.5,beta:3.0,p:-0.7,q:3.5,zcut:2.0}
     ldat1={phi0:1.0,lo:1.0,alpha:1.0,beta:1.0,p:0.0,q:0.0,zcut:1.0}
     ldat2={phi0:-5.0,lo:9.0,alpha:0.0,beta:0.0,p:-6.0,q:2.0,zcut:0.0}
     ldat3={phi0:5.0,lo:11.0,alpha:2.0,beta:5.0,p:2.0,q:8.0,zcut:0.0}
     ldat4={phi0:0.0,lo:0.0,alpha:0.0,beta:0.0,p:0.3,q:0.1,zcut:0.0}
     ldata=[ldat0,ldat1,ldat2,ldat3,ldat4]
                                ;the fixed values: =1 if held fixed, =0 if variable
     
     cdat = {a0:0.0,fixed:0.0,amin:0.0,amax:1.0,sigma:0.1}
     sdat = {area:10.0,zmin:0.0,zmax:5.0,dz:0.1,runs:1.e3}
     
     CD, Current=thisdir
     files = {ofile:thisdir+'/observation.save', $
              mfile:thisdir+'/model.fits', $
              sedfile:thisdir+'/sf_templates.fits', $
              oname:thisdir+'/output.fits' }
     
     msettings = {$
                 nchain:5,$
                 tmax:20.0,$
                 acceptpct:0.25, $
                 pct_range:0.05, $
                 conv_conf:0.05, $
                 conv_rmax:1.05, $
                 conv_step:20, $
                 burn_step:10, $
                 burn_ratio:10}
  endelse
  
  bname = ["Band 1","Band 2","Band 3"]
  ocols = ["Wavelength (um)","Flux limit (mJy)","Standard Error (mJy)"]
  f = ['(e9.2)','(f7.4)','(f7.4)']
  fmt = [[f],[f],[f]]
  
;The Survey properties table
  lo = widget_label(info.obs_table,value="Survey Properties")
  info.ot = widget_table(info.obs_table,value=bands,column_labels=ocols,row_labels=bname,uvalue='ot',/editable,alignment=1,column_widths=[100,100,150],format=fmt,scr_xsize=425,scr_ysize=95)
  
  f2a = ['(f5.2)','(f5.2)','(f5.2)','(f5.2)','(f5.2)','(f5.2)','(f5.2)']
  f2b = ['(i)','(i)','(i)','(i)','(i)','(i)','(i)']
  fmt2 = [[f2a],[f2b],[f2a],[f2a],[f2a]]
  fmt3 = ['(f5.2)','(i)','(f5.2)','(f5.2)','(f5.2)']
  lrows=["Initial","Fixed","Min","Max","Sigma"]
  
;Luminosity Function Parameters
  l1 = widget_label(info.lum_table,value="Luminosity Function Parameters")
  info.t1 = widget_table(info.lum_table,value=ldata,column_labels=tag_names(ldata),row_labels=lrows,uvalue='t1',/editable,alignment=1,format=fmt2,scr_xsize=537,scr_ysize=132)
  
  l3 = widget_label(info.sed_table,value="SED Evolution Parameters")
  info.t3 = widget_table(info.sed_table,value=cdat,column_labels=lrows,row_labels=["Color Exp"],uvalue='t3',/editable,alignment=1,format=fmt3,scr_xsize=404,scr_ysize=55)
  
  tcols = ["Area (sdeg)","Z Min","Z Max","Z Binsize","Run Number"]
;Simulation Parameters
  l2 = widget_label(info.sim_table,value="Simulation Settings")
  info.t2 = widget_table(info.sim_table,value=sdat,column_labels=tcols,/no_row_headers,uvalue='t2',/editable,alignment=1,format=['(f5.2)','(f5.2)','(f5.2)','(f5.2)','(e9.2)'],column_widths=[100,100,100,100,100],scr_xsize=505,scr_ysize=55)
  
  
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
                           ysize=info.magnification*info.ysize1 $
                           ,uvalue="DRAW_WINDOW3",retain=2 $
                           ,/button_events, keyboard_events=1,/tracking_events)
  
;initialize widget and establish control
  widget_control,/realize,info.base,xoffset=0,yoffset=0,Set_UValue=theObject
  widget_control, info.draw3, get_value=win_id3 & info.win_id3=win_id3
  xmanager,'simulation',info.base,/NO_BLOCK
  wset,info.win_id3
  
;plot SEDs (generalize to any SED template file)
  set_plot,'x'
  device,decomposed=0
  plot_settings,plot_type='x'
  if file_test(files.sedfile) then begin
     templ=mrdfits(files.sedfile,/silent) 
     loadct,1,/silent
     plot,templ[*,0],templ[*,1],/xlog,/ylog,yrange=[1.d20,1.d28],ystyle=1,xtitle=TeXtoIDL('\lambda [\mum]'),ytitle=TeXtoIDL('L_{\nu} [W/Hz]')
     for ipl=1,13 do oplot,templ[*,0],templ[*,ipl+1]
  endif else begin
     print,'File '+files.sedfile+' does not exist, skipping plot'
  endelse
  
END

;widget event handling routine
PRO simulation_event,ev
  COMMON simulation_com,info,ldata,ldat0,sdat,cdat,bands,msettings,files

  ; get event identifier
  widget_control,ev.id,get_uvalue=uvalue

  CASE uvalue OF
     'save' : save,ldata,ldat0,bands,sdat,cdat,msettings,files,filename='params.save' ;save parameters
     'go'  : begin              ;save settings, intialize FITS file, and pass to C++ fitting routine
        save,ldata,ldat0,bands,sdat,cdat,msettings,files,filename='params.save' ;save parameters
        
                                ;make observation FITS file
        widget_control,info.ot,get_value=bvals
        wave = [bvals[0].wave,bvals[1].wave,bvals[2].wave]
        flux_min = [bvals[0].fmin,bvals[1].fmin,bvals[2].fmin]
        flux_err = [bvals[0].ferr,bvals[1].ferr,bvals[2].ferr]

        if (file_test(files.ofile)) then begin
           print,files.ofile
           restore,files.ofile
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

        widget_control,info.t1,get_value=lparam
        pars = lparam(0)
        fixed = lparam(1)
        min = lparam(2)
        max = lparam(3)
        sigma = lparam(4)
        widget_control,info.t2,get_value=sparam
        sdat = sparam
        widget_control,info.t3,get_value=cparam
        cdat = cparam

        sxaddpar,hdr2,'DATE',systime(),'Date of creation'
        sxaddpar,hdr2,'PHI0',pars.phi0,'Luminosity Function Normalization'
        sxaddpar,hdr2,'PHI0_FIX',fixed.phi0,'Fix Phi0 (Y=1/N=0)'
        sxaddpar,hdr2,'PHI0_MIN',min.phi0,'Minimum Phi0 value'
        sxaddpar,hdr2,'PHI0_MAX',max.phi0,'Maximum Phi0 value'
        sxaddpar,hdr2,'PHI0_DP',sigma.phi0,'Sigma Phi0 value'
        sxaddpar,hdr2,'L0',pars.lo,'Luminosity Function Knee'
        sxaddpar,hdr2,'L0_FIX',fixed.lo,'Fix L0 (Y=1/N=0)'
        sxaddpar,hdr2,'L0_MIN',min.lo,'Minimum L0 value'
        sxaddpar,hdr2,'L0_MAX',max.lo,'Maximum L0 value'
        sxaddpar,hdr2,'L0_DP',sigma.lo,'Sigma L0 value'
        sxaddpar,hdr2,'ALPHA',pars.alpha,'Luminosity Function upper slope'
        sxaddpar,hdr2,'ALPHA_FIX',fixed.alpha,'Fix Alpha (Y=1/N=0)'
        sxaddpar,hdr2,'ALPHA_MIN',min.alpha,'Minimum Alpha value'
        sxaddpar,hdr2,'ALPHA_MAX',max.alpha,'Maximum Alpha value'
        sxaddpar,hdr2,'ALPHA_DP',sigma.alpha,'Sigma Alpha value'
        sxaddpar,hdr2,'BETA',pars.beta,'Luminosity Function lower slope'
        sxaddpar,hdr2,'BETA_FIX',fixed.beta,'Fix Beta (Y=1/N=0)'
        sxaddpar,hdr2,'BETA_MIN',min.beta,'Minimum Beta value'
        sxaddpar,hdr2,'BETA_MAX',max.beta,'Maximum Beta value'
        sxaddpar,hdr2,'BETA_DP',sigma.beta,'Sigma Beta value'
        sxaddpar,hdr2,'P',pars.p,'Luminosity Function PHI evolution term'
        sxaddpar,hdr2,'P_FIX',fixed.p,'Fix P (Y=1/N=0)'
        sxaddpar,hdr2,'P_MIN',min.p,'Minimum P value'
        sxaddpar,hdr2,'P_MAX',max.p,'Maximum P value'
        sxaddpar,hdr2,'P_DP',sigma.p,'Sigma P value'
        sxaddpar,hdr2,'Q',pars.q,'Luminosity Function L evolution term'
        sxaddpar,hdr2,'Q_FIX',fixed.q,'Fix Q (Y=1/N=0)'
        sxaddpar,hdr2,'Q_MIN',min.q,'Minimum Q value'
        sxaddpar,hdr2,'Q_MAX',max.q,'Maximum Q value'
        sxaddpar,hdr2,'Q_DP',sigma.q,'Sigma Q value'
        sxaddpar,hdr2,'ZCUT',pars.zcut,'Luminosity Function z evolution limit'
        sxaddpar,hdr2,'ZCUT_FIX',fixed.zcut,'Fix ZCUT (Y=1/N=0)'
        sxaddpar,hdr2,'ZCUT_MIN',min.zcut,'Minimum ZCUT value'
        sxaddpar,hdr2,'ZCUT_MAX',max.zcut,'Maximum ZCUT value'
        sxaddpar,hdr2,'ZCUT_DP',sigma.zcut,'Sigma ZCUT value'
        sxaddpar,hdr2,'CEXP',cdat.a0,'Intrinsic luminosity evolution term'
        sxaddpar,hdr2,'CEXP_FIX',cdat.fixed,'Fix CEXP (Y=1/N=0)'
        sxaddpar,hdr2,'CEXP_MIN',cdat.amin,'Minimum CEXP value'
        sxaddpar,hdr2,'CEXP_MAX',cdat.amax,'Maximum CEXP value'
        sxaddpar,hdr2,'CEXP_DP',cdat.sigma,'Sigma CEXP value'
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
        args = 'observation.fits '+files.mfile+' '+files.sedfile+' '+files.oname
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
     'info': fit_info
     ELSE:
  ENDCASE

END

pro settings
  
  COMMON simulation_com,info,ldata,ldat0,sdat,cdat,bands,msettings,files
  
  smain = widget_base(title="Simulation Settings",/column)
  body = widget_base(smain,/row,/align_center)
  sbase = widget_base(body,/column,/align_center)
  fbase = widget_base(body,/column,/align_center)

  l1 = widget_label(sbase,value="Monte Carlo Runtime Parameters")
  info.tset = widget_table(sbase,value=[msettings],row_labels=tag_names(msettings),column_labels=["Value"],/editable,alignment=1,/column_major,uvalue='settings')

  info.dprint = widget_droplist(fbase,value=["yes","no"],title="Enable Verbose Simulation Output",uvalue="print")
  widget_control,info.dprint,set_droplist_select=info.print
  
  info.obsname = fsc_fileselect(fbase,/NoMaxSize,LabelName='Observation Save File')
  widget_control,info.obsname,set_value=files.ofile
  
  info.sfile = fsc_fileselect(fbase,ObjectRef=sedObject,/NoMaxSize,LabelName='SED Templates File')
  widget_control,info.sfile,set_value=files.sedfile

  info.mname = fsc_fileselect(fbase,/NoMaxSize,LabelName='Model Fits File')
  widget_control,info.mname,set_value=files.mfile

  info.oname = fsc_fileselect(fbase,/NoMaxSize,LabelName='Output Fits File')
  widget_control,info.oname,set_value=files.oname

  btns = widget_base(smain,/row,/align_left)
  save_btn = widget_button(btns,uvalue='save',Value="Save and Close")
  cancel_btn = widget_button(btns,uvalue='cancel',Value="Cancel")

  widget_control,smain,/realize
  xmanager,'settings',smain,/no_block

end

pro settings_event,ev
  
  COMMON simulation_com,info,ldata,ldat0,sdat,cdat,bands,msettings,files
  
  widget_control,ev.id,get_uvalue=uvalue

  switch uvalue of
     'save' : begin
        widget_control,info.tset,get_value=msettings
        info.print = widget_info(info.dprint,/droplist_select)
        widget_control,info.sfile,get_value=temp
        files.sedfile=temp
        widget_control,info.obsname,get_value=temp
        files.ofile=temp
        widget_control,info.mname,get_value=temp
        files.mfile=temp
        widget_control,info.oname,get_value=temp
        files.oname=temp
     end
     'cancel': begin
        widget_control,ev.top,/destroy
     end
  endswitch

end

pro fit_info
  
  main = widget_base(title="Information about Fitting Methodology",/column)
  t1 = widget_text(main,/wrap,value=['Luminosity Function Fitter Version 1.0' $
  ,'Noah Kurinsky and Anna Sajina' $
  ,'Department of Physics and Astronomy, Tufts University, Medford MA' $
  ,'Last Updated 5/1/13' $
  ,'' $
  ,'Notes on using this interface:' $
  ,'- To allow a parameter to be fitted, set the value of the "fixed" row for that parameter to "0"' $
  ,'- The ranges for the luminosity function parameters are only important for non-fixed parameters' $
  ,'- The "sigma" row refers to the proposed distribution width (stdev of gaussian monte carlo steps'],xsize=80,ysize=20)
  bthold = widget_base(main,/row,/align_center)
  bt = widget_button(bthold,uvalue='close',value='Close',xsize=50,ysize=25)

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
  
  COMMON simulation_com,info,ldata,ldat0,sdat,cdat,bands,msettings,files

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
  c1mean = []
  c2mean = []
  c3mean = []
  c1plus = []
  c2plus = []
  c3plus = []
  c1minus = []
  c2minus = []
  c3minus = []
  c1size = n_elements(count_dists[0].dnds250)
  c2size = n_elements(count_dists[0].dnds350)
  c3size = n_elements(count_dists[0].dnds500)

  for i=0,c1size-1 do begin
     dnds = c1[i,*]
     dnds = dnds[sort(dnds)]
     dnds = dnds[where(dnds gt 0)]
     pi = plusfrac*n_elements(dnds)
     mi = minusfrac*n_elements(dnds)
     c1mean = [c1mean,mean(dnds)]
     c1plus = [c1plus,dnds[pi]]
     c1minus = [c1minus,dnds[mi]]
  endfor

  for i=0,c2size-1 do begin
     dnds = c2[i,*]
     dnds = dnds[sort(dnds)]
     dnds = dnds[where(dnds gt 0)]
     pi = plusfrac*n_elements(dnds)
     mi = minusfrac*n_elements(dnds)
     c2mean = [c2mean,mean(dnds)]
     c2plus = [c2plus,dnds[pi]]
     c2minus = [c2minus,dnds[mi]]
  endfor

  for i=0,c3size-1 do begin
     dnds = c3[i,*]
     dnds = dnds[sort(dnds)]
     dnds = dnds[where(dnds gt 0)]
     pi = plusfrac*n_elements(dnds)
     mi = minusfrac*n_elements(dnds)
     c3mean = [c3mean,mean(dnds)]
     c3plus = [c3plus,dnds[pi]]
     c3minus = [c3minus,dnds[mi]]
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
  
  chis = []
  cpts = where(strmatch(alltags,'CHISQ*',/FOLD_CASE) eq 1)
  for ci=0,n_elements(cpts)-1 do begin
     chis = [chis,res.(cpts[ci])]
  endfor

  accept = []
  apts = where(strmatch(alltags,'ACPT*',/FOLD_CASE) eq 1)
  for ai=0,n_elements(apts)-1 do begin
     accept = [accept,res.(apts[ai])]
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
     p = []
     ppts = where(strmatch(alltags,tags[x]+'*',/FOLD_CASE) eq 1)
     for pi=0,n_elements(ppts)-1 do begin
        p = [p,res.(ppts[pi])]
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

           q = []
           qpts = where(strmatch(alltags,tags[y]+'*',/FOLD_CASE) eq 1)
           for qi=0,n_elements(qpts)-1 do begin
              q = [q,res.(qpts[qi])]
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
                 chi_plot = mean(chis[gpts])
                 
                 if(gpts[0] ne -1) then begin
                    polyfill,xfill,yfill,color=color_scale*(chi_max - alog(chi_plot))
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

  COMMON simulation_com,info,ldata,ldat0,sdat,cdat,bands,msettings,files

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
  c1mean = []
  c2mean = []
  c3mean = []
  c1plus = []
  c2plus = []
  c3plus = []
  c1minus = []
  c2minus = []
  c3minus = []
  c1size = n_elements(count_dists[0].dnds250)
  c2size = n_elements(count_dists[0].dnds350)
  c3size = n_elements(count_dists[0].dnds500)

  for i=0,c1size-1 do begin
     dnds = c1[i,*]
     dnds = dnds[sort(dnds)]
     dnds = dnds[where(dnds gt 0)]
     pi = plusfrac*n_elements(dnds)
     mi = minusfrac*n_elements(dnds)
     c1mean = [c1mean,mean(dnds)]
     c1plus = [c1plus,dnds[pi]]
     c1minus = [c1minus,dnds[mi]]
  endfor

  for i=0,c2size-1 do begin
     dnds = c2[i,*]
     dnds = dnds[sort(dnds)]
     dnds = dnds[where(dnds gt 0)]
     pi = plusfrac*n_elements(dnds)
     mi = minusfrac*n_elements(dnds)
     c2mean = [c2mean,mean(dnds)]
     c2plus = [c2plus,dnds[pi]]
     c2minus = [c2minus,dnds[mi]]
  endfor

  for i=0,c3size-1 do begin
     dnds = c3[i,*]
     dnds = dnds[sort(dnds)]
     dnds = dnds[where(dnds gt 0)]
     pi = plusfrac*n_elements(dnds)
     mi = minusfrac*n_elements(dnds)
     c3mean = [c3mean,mean(dnds)]
     c3plus = [c3plus,dnds[pi]]
     c3minus = [c3minus,dnds[mi]]
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
  resmod[neg] = 0
  resobs[pos] = 0

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
  COMMON simulation_com,info,ldata,ldat0,sdat,cdat,bands,msettings,files

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

  COMMON simulation_com,info,ldata,ldat0,sdat,cdat,bands,msettings,files

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
  
  chis = []
  cpts = where(strmatch(alltags,'CHISQ*',/FOLD_CASE) eq 1)
  for ci=0,n_elements(cpts)-1 do begin
     chis = [chis,res.(cpts[ci])]
  endfor

  accept = []
  apts = where(strmatch(alltags,'ACPT*',/FOLD_CASE) eq 1)
  for ai=0,n_elements(apts)-1 do begin
     accept = [accept,res.(apts[ai])]
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
     p = []
     ppts = where(strmatch(alltags,tags[x]+'*',/FOLD_CASE) eq 1)
     for pi=0,n_elements(ppts)-1 do begin
        p = [p,res.(ppts[pi])]
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

           q = []
           qpts = where(strmatch(alltags,tags[y]+'*',/FOLD_CASE) eq 1)
           for qi=0,n_elements(qpts)-1 do begin
              q = [q,res.(qpts[qi])]
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
                 chi_plot = mean(chis[gpts])
                 
                 if(gpts[0] ne -1) then begin
                    polyfill,xfill,yfill,color=color_scale*(chi_max - alog(chi_plot))
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
  COMMON simulation_com,info,ldata,ldat0,sdat,cdat,bands,msettings,files

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

pro plot_settings,plot_type=plot_type
  cleanplot,/silent
  
  if not keyword_set(plot_type) then plot_type='reset'
  
  if strlowcase(plot_type) eq 'ps' then begin
     !p.thick=5
     !x.thick=5
     !y.thick=5
     !p.charthick=5
     !p.charsize=1.5
     cleanplot,/showonly
  endif else if strlowcase(plot_type) eq 'x' then begin
     !p.background=255
     !p.color=0
     cleanplot,/showonly
  endif else begin
     print,"Returning to default plot settings"
  endelse
end