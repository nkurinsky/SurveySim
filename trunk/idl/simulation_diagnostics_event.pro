pro simulation_diagnostics_event,ev
  COMMON simulation_com

  widget_control,ev.id,get_uvalue=uvalue

  case uvalue of
     'graphs': begin
        widget_control,ev.top,/destroy
        simulation_results
     end
     'refresh': begin
        widget_control,ev.top,get_uvalue=mnum
        widget_control,ev.top,/destroy
        simulation_diagnostics
     end
     'close': widget_control,ev.top,/destroy
     'quit': begin
        widget_control,ev.top,/destroy
        widget_control,info.base,/destroy
     end
  endcase

end
