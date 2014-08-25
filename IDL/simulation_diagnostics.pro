pro simulation_diagnostics

  COMMON simulation_com

  if(n_elements(parameters) eq 0) then begin
     if(keyword_set(savefile)) then begin
        temp = load_parameters(savefile)
     endif else begin
        temp = load_parameters()
     endelse
     parameters = temp
  endif
  
  filename = parameters.files.oname

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
  xmanager,'simulation_diagnostics',dmain,/no_block
  
  set_plot,'x'
  plot_settings,plot_type='x'
  device,decomposed=0
  loadct,0,/silent
  
; chain read operations 
  widget_control,resplot,get_value=index
  wset,index
  plot_fit_results,filename

  widget_control,chisqr,get_value=index
  wset,index
  plot_chisq_trend,filename
  
  widget_control,chidist,get_value=index
  wset,index
  plot_chisq_hist,filename

  widget_control,conv,get_value=index
  wset,index
  plot_convergence,filename

end
