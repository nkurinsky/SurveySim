pro ask_for_file_event,ev

  common ask_for_file_com,new_filename,dialog

  widget_control,ev.id,get_uvalue=uvalue

  switch uvalue of
     'select' : begin
        widget_control,dialog,get_value=new_filename
     end
     'skip': begin
        widget_control,ev.top,/destroy
     end
  endswitch
  
end
