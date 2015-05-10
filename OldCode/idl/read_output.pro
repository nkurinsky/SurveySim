pro read_output,savefile=savefile
  
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
  sdat = parameters.surveyData

  print,filename
    
  set_plot,'ps'
  plot_settings,plot_type='ps'

  device,filename='redshift_dist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  plot_redshift,filename,sdat.zmin,sdat.zmax,sdat.dz
  device,/close

  device,filename='band1_counts.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  plot_counts,filename,1
  device,/close
  
  device,filename='band2_counts.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  plot_counts,filename,2
  device,/close

  device,filename='band3_counts.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  plot_counts,filename,3
  device,/close

  device,filename='comp_color_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color  
  plot_simulation_hist,filename,0
  device,/close

  device,filename='model_color_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  plot_simulation_hist,filename,1
  device,/close

  device,filename='obs_color_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  plot_simulation_hist,filename,2
  device,/close

  device,filename='sim_lumfunct.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated 
  plot_lumfunct
  device,/close

; chain read operations 

  set_plot,'ps'
  plot_settings,plot_type='ps'
  device,filename='fit_results.eps',xsize=10,ysize=8,/inches,/times,set_font='Times-Roman',/color,/encapsulated
  plot_fit_results,filename
  device,/close
  
  device,filename='chisq_v_run.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated
  plot_chisq_trend,filename
  device,/close

  device,filename='chisq_hist.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated
  plot_chisq_hist,filename
  device,/close

  device,filename='convergence.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated
  plot_convergence,filename
  device,/close

  if(file_test('plots')) then begin
     file_delete,'plots',/recursive
  endif
  file_mkdir,'plots'
  file_move,'*.eps','plots'

end
