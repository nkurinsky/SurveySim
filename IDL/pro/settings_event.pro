pro settings_event,ev
  
  COMMON simulation_com
  
  widget_control,ev.id,get_uvalue=uvalue

  switch uvalue of
     'save' : begin
        widget_control,info.tset,get_value=msettings
        info.print = widget_info(info.dprint,/droplist_select)
        
        widget_control,info.sfile,get_value=files.sedfile
        widget_control,info.obsname,get_value=files.ofile
        widget_control,info.mname,get_value=files.mname
        widget_control,info.oname,get_value=files.oname
     end
     'cancel': begin
        widget_control,ev.top,/destroy
     end
  endswitch

end
