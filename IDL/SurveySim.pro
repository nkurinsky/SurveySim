;==================================================================
; Writen by Noah Kurinsky, version recent as of 8/21/14
;==================================================================

COMMON simulation_com,info,parameters

@SurveySimIncludes
@SurveySim_event

pro SurveySim,savefile=savefile

  COMMON simulation_com,info,parameters
  
  spawn,'pwd',thisdir
  plot_settings                 ;will reset the global settings
  
  ;==================================================================
  ;the INFO structure holds the key widget control parameters as well as 
  ;the basic simulation settings
  ;------------------------------------------------------------------
  if(keyword_set(savefile)) then begin
     temp=load_parameters(savefile)
  endif else begin
     temp=load_parameters()
  endelse
  parameters = temp
  info = make_info_struct(n_elements(parameters.lumpars))

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
  filter_names = filter_list("/usr/local/surveysim/filters/filterlib.txt")
  filter_names = ["Select Filter",filter_names]
  
  ;The filter properties table
  lo = widget_label(info.obs_table,value="Survey Properties")
  obs_row = widget_base(info.obs_table,/row)
  filter_dialogs = widget_base(obs_row,/column,/base_align_bottom,/align_bottom)
  info.fd1 = widget_combobox(filter_dialogs, $
                             value=filter_names,$
                             uvalue="fd1")
  info.fd2 = widget_combobox(filter_dialogs, $
                             value=filter_names,$
                             uvalue="fd2")
  info.fd3 = widget_combobox(filter_dialogs, $
                             value=filter_names,$
                             uvalue="fd3")
  info.ot = widget_table(obs_row,$
                         uvalue='ot',/editable,alignment=1,$
                         value=parameters.filters.properties,$
                         column_labels=["Flux limit (mJy)","Standard Error (mJy)"],$
                         /no_row_headers,$
                         column_widths=150, row_heights=28,$
                         format='(f5.2)',$
                         scr_xsize=304,scr_ysize=120)

  widget_control, info.fd1, set_combobox_select=parameters.filters[0].filter_id
  widget_control, info.fd2, set_combobox_select=parameters.filters[1].filter_id
  widget_control, info.fd3, set_combobox_select=parameters.filters[2].filter_id

  ;Survey Model Parameters
  l1 = widget_label(info.lum_table,value="Survey Model Parameters")
  lumrow = widget_base(info.lum_table,/row,/align_center)
  info.t1 = widget_table(lumrow, $
                         value=parameters.lumpars.pars,$
                         row_labels=parameters.lumpars.name,$
                         column_labels=tag_names(parameters.lumpars.pars),$
                         column_widths=100, row_heights=28, $
                         uvalue='t1',$
                         /editable,alignment=1,$
                         format='(f5.2)',$
                         scr_xsize=374,scr_ysize=36+n_elements(parameters.lumpars.pars)*28)
  
  fixvalues = parameters.lumpars.fixed
  fnum = n_elements(fixvalues)
  fixcol = widget_base(lumrow,/column,/align_bottom,/base_align_bottom)
  for i=0,fnum-1 do begin
     info.fixinfo[i] = widget_combobox(fixcol, $
                                       value=["Fitted","Fixed"],$
                                       uvalue="fix"+strtrim(string(i),1))
     widget_control, info.fixinfo[i], set_combobox_select=fixvalues[i]
  endfor
  
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
  widget_control, info.draw3, get_value= temp
  info.win_id3 = temp

  xmanager,'SurveySim',info.base,/NO_BLOCK
  wset,info.win_id3

  ;get new sedfile if file does not exist (option to ignore too)
  if not file_test(parameters.files.sedfile) then begin
     parameters.files.sedfile = ask_for_file("SED File "+ $
                                            parameters.files.sedfile+ $
                                            " not found, please enter new SED file", $
                                            parameters.files.sedfile)
  endif
  
  plot_seds
  
END
