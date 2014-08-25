pro settings
  
  COMMON simulation_com
  
  smain = widget_base(title="Simulation Settings",/column)
  body = widget_base(smain,/row,/align_center)
  sbase = widget_base(body,/column,/align_center)
  fbase = widget_base(body,/column,/align_center)

  l1 = widget_label(sbase,value="MCMC Runtime Parameters")
  info.tset = widget_table(sbase,value=[parameters.msettings],row_labels=tag_names(parameters.msettings),column_labels=["Value"],/editable,alignment=1,/column_major,uvalue='settings',scr_xsize=70,scr_ysize=175,/no_column_headers)

  info.dprint = widget_droplist(fbase,value=["yes","no"],title="Enable Verbose Simulation Output",uvalue="print")
  widget_control,info.dprint,set_droplist_select=parameters.print
  
  info.obsname = fsc_fileselect(fbase,/NoMaxSize,$
                                LabelName='Observation Save File',$
                                SelectTitle='Select Observation FITS File...')
  widget_control,info.obsname,set_value=parameters.files.ofile
  
  info.sfile = fsc_fileselect(fbase,/NoMaxSize,$
                              LabelName='SED Templates File',$
                              SelectTitle='Select SED FITS File...')
  widget_control,info.sfile,set_value=parameters.files.sedfile

  info.mname = fsc_fileselect(fbase,/NoMaxSize,$
                              LabelName='Model Fits File',$
                              SelectTitle='Select Model FITS File...')
  widget_control,info.mname,set_value=parameters.files.mfile

  info.oname = fsc_fileselect(fbase,/NoMaxSize,$
                              LabelName='Output Fits File',$
                              SelectTitle='Select Output FITS File...')
  widget_control,info.oname,set_value=parameters.files.oname

  btns = widget_base(smain,/row,/align_left)
  save_btn = widget_button(btns,uvalue='save',Value="Save and Close")
  cancel_btn = widget_button(btns,uvalue='cancel',Value="Cancel")
  
  widget_control,smain,/realize
  xmanager,'settings',smain

end
