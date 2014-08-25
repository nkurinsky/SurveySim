pro settings_event,ev
  
  COMMON simulation_com
  
  widget_control,ev.id,get_uvalue=uvalue

  switch uvalue of
     'save' : begin
        widget_control,info.tset,get_value=temp
        parameters.msettings = temp[0]
        parameters.print = widget_info(info.dprint,/droplist_select)
        
        widget_control,info.sfile,get_value=temp
        parameters.files.sedfile = temp
        widget_control,info.obsname,get_value=temp
        parameters.files.ofile=temp
        widget_control,info.mname,get_value=temp
        parameters.files.mfile=temp
        widget_control,info.oname,get_value=temp
        parameters.files.oname=temp
     end
     'cancel': begin
        widget_control,ev.top,/destroy
     end
  endswitch

end
