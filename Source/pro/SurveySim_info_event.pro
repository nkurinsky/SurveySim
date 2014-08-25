pro SurveySim_info_event,ev

  widget_control,ev.id,get_uvalue=uvalue

  case uvalue of
     'close': begin
        widget_control,ev.top,/destroy
     end
  endcase

end
