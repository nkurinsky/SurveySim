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
  if not file_test(parameters.files.sedfile) then begin
     parameters.files.sedfile = ask_for_file("SED File "+ $
                                            parameters.files.sedfile+ $
                                            " not found, please enter new SED file", $
                                            parameters.files.sedfile)
  endif
  
  ;plot SEDs
  set_plot,'x'
  device,decomposed=0
  plot_settings,plot_type='x'
  if file_test(parameters.files.sedfile) then begin

     templ=mrdfits(parameters.files.sedfile,/silent) 
     loadct,1,/silent

     plot,templ[*,0],templ[*,1],/xlog,/ylog,$
          yrange=[1.d20,1.d28],ystyle=1,ytitle=TeXtoIDL('L_{\nu} [W/Hz]'),$
          xtitle=TeXtoIDL('\lambda [\mum]')

     for ipl=1,13 do oplot,templ[*,0],templ[*,ipl+1]
  endif else begin
     print,'File '+parameters.files.sedfile+' does not exist, skipping plot'
  endelse
  
END

;widget event handling routine
PRO SurveySim_event,ev
  COMMON simulation_com
  
  ; get event identifier
  widget_control,ev.id,get_uvalue=uvalue

  CASE uvalue OF
     'save' : save,parameters,filename='params.save' ;save parameters
     'go'  : begin              ;save settings, intialize FITS file, and pass to C++ fitting routine
        save,parameters,filename='params.save' ;save parameters
        
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












