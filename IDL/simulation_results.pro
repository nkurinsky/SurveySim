pro simulation_results,savefile=savefile

  COMMON simulation_com

  if(n_elements(parameters) eq 0) then begin
     if(keyword_set(savefile)) then begin
        temp = load_parameters(savefile)
     endif else begin
        temp = load_parameters()
     endelse
     parameters = temp
  endif

  size_screen=get_screen_size()
  size_screen_alt = size_screen*0.85
  size_screen = size_screen*0.8
  gmain = widget_base(title='Simulation Output',/column,xsize=size_screen[0],ysize=size_screen_alt[1])
  r1 = widget_base(gmain,/row)
  r2 = widget_base(gmain,/row)
  r2b = widget_base(gmain,/row)
  r3 = widget_base(gmain,/row)

  widget_control,gmain,set_uvalue=mnum
  xdim = fix(size_screen[0]/3.0)
  ydim = fix(size_screen[1]/3.0)

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
  xmanager,'simulation_results',gmain,/no_block

  filename = parameters.files.oname
  sdat = parameters.surveyData
  
  set_plot,'x'
  plot_settings,plot_type='x'
  loadct,0,/silent
  device,decomposed=0

  widget_control,lumfunct,get_value=index
  wset,index
  plot_lumfunct

  widget_control,redshift,get_value=index
  wset,index
  plot_redshift,filename,sdat.zmin,sdat.zmax,sdat.dz

  widget_control,models,get_value=index
  wset,index
  plot_lumdist,filename

  widget_control,dcount1,get_value=index
  wset,index
  plot_counts,filename,1

  widget_control,dcount2,get_value=index
  wset,index
  plot_counts,filename,2

  widget_control,dcount3,get_value=index
  wset,index
  plot_counts,filename,3

  widget_control,sim_colors,get_value=index
  wset,index
  plot_simulation_hist,filename,1

  widget_control,obs_colors,get_value=index
  wset,index
  plot_simulation_hist,filename,2

  widget_control,comp_colors,get_value=index
  wset,index
  plot_simulation_hist,filename,0

end
