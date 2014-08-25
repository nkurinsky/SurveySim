;==================================================================
; Writen by Noah Kurinsky, version recent as of 8/21/14
;==================================================================

PRO Run_Simulation
  COMMON simulation_com

  save,parameters,filename='params.save'
  
  filter_names = filter_list("/usr/local/surveysim/filters/filterlib.txt")
  flux_min = parameters.filters.properties.min
  flux_err = parameters.filters.properties.err
  filterids = parameters.filters.filter_id

  obs_fitsfile = 'observation.fits'
  spawn,'cp '+parameters.files.ofile+' '+obs_fitsfile
  
  hdr = headfits(obs_fitsfile)
  hnum = sxpar(hdr,"FHDU")
  hdr = headfits(obs_fitsfile,exten=hnum)

  for i=0,2 do begin
     str_ind = strtrim(string(i),1)
     if(filterids[i] ne 0) then begin
        sxaddpar, hdr, 'F'+str_ind+'MIN',flux_min[i],'Flux cutoff, column '+str_ind
        sxaddpar, hdr, 'F'+str_ind+'ERR',flux_err[i],'Flux error, column '+str_ind
        sxaddpar, hdr, 'F'+str_ind+'FILT',filter_names[filterids[i]],'Filter name, column '+str_ind
     endif else begin
        print,"Please select filter for column "+str_ind
        return
     endelse
  endif
     
  modfits,obs_fitsfile,0,hdr,exten_no=hnum
  
  ;write model keywords
  sxaddpar,hdr2,'DATE',systime(),'Date of creation'
  
  for i=0,n_elements(parameters.lumpars)-1 do begin
     strnum = strtrim(string(i),1)
     parname=parameters.lumpars[i].name
     sxaddpar,hdr2,parname,parameters.lumpars[i].pars.value,'Luminosity Function Parameter '+strnum
     sxaddpar,hdr2,parname+'_FIX',parameters.lumpars[i].fixed,'Fix '+parname+' (Y=1/N=0)'
     sxaddpar,hdr2,parname+'_MIN',parameters.lumpars[i].min,'Minimum '+parname+' value'
     sxaddpar,hdr2,parname+'_MAX',parameters.lumpars[i].max,'Maximum '+parname+' value'
  endfor
  
  sxaddpar,hdr2,'RUNS',parameters.surveyData.runs,'Number of Runs'
  sxaddpar,hdr2,'ZMIN',parameters.surveyData.zmin,'Minimum Redshift Value'
  sxaddpar,hdr2,'ZMAX',parameters.surveyData.zmax,'Maximum Redshift Value'
  sxaddpar,hdr2,'DZ',parameters.surveyData.dz,'Redshit Bin Width'
  sxaddpar,hdr2,'AREA',parameters.surveyData.area,'Observed Solid Angle'
  
  sxaddpar,hdr2,'NCHAIN',parameters.msettings.nchain,'Chain Number'
  sxaddpar,hdr2,'TMAX',parameters.msettings.tmax,'Starting Anneal Temperature'
  sxaddpar,hdr2,'ANN_PCT',parameters.msettings.acceptpct,'Ideal Acceptance Percentage'
  sxaddpar,hdr2,'CONV_CONF',parameters.msettings.conv_conf,'Convergence CI Setting'
  sxaddpar,hdr2,'CONV_RMAX',parameters.msettings.conv_rmax,'Convergence Rmax Criterion'
  sxaddpar,hdr2,'CONV_STEP',parameters.msettings.conv_step,'Iterations between convergence checks'
  sxaddpar,hdr2,'BURN_STEP',parameters.msettings.burn_step,'Iterations between anneal calls in burn-in'
  sxaddpar,hdr2,'BURNVRUN',parameters.msettings.burn_ratio,'Ratio of normal to burn-in steps'
  sxaddpar,hdr2,'ANN_RNG',parameters.msettings.pct_range,'Range within which to maintain acceptance'
  sxaddpar,hdr2,'PRINT',parameters.print,'Whether to Print Debug MSGs'
  
  mwrfits,[0],parameters.files.mfile,hdr2,/create
  
                                ;Run the actual simulation
  command = 'fitter '+ $
            obs_fitsfile+' '+ $
            parameters.files.mfile+' '+ $
            parameters.files.sedfile+' '+$
            parameters.files.oname
  print,'Running Command "'+command+'"'
  spawn,command,exit_status=result
  
  if(result eq 0) then begin
     read_output,parameters.files.oname
     simulation_results,parameters.files.oname
  endif else begin
     print,"Fitting Failure"
  endelse
END

PRO SurveySim_event,ev
  COMMON simulation_com
  
  ; get event identifier
  widget_control,ev.id,get_uvalue=uvalue
  case_value = strmid(uvalue,0,3)

  CASE case_value OF
     'sav' : save,parameters,filename='params.save' ;save parameters
     'go'  : Run_Simulation              
     'set': begin
        settings
        plot_seds
     end
     'dia': simulation_diagnostics
     'rep': begin
        read_output,parameters.files.oname
        simulation_results,parameters.files.oname
     end
     'qui': widget_control,ev.top,/destroy
     'ot': widget_control,info.ot,get_value=parameters.filters.properties
     'fd1': widget_control,info.fd1, get_combobox_select=parameters.filters[0].filter_id
     'fd2': widget_control,info.fd2, get_combobox_select=parameters.filters[1].filter_id
     'fd3': widget_control,info.fd3, get_combobox_select=parameters.filters[2].filter_id
     't1': widget_control,info.t1,get_value=parameters.lumpars.pars
     't2': widget_control,info.t2,get_value=parameters.surveyData
     'inf': SurveySim_info
     'fix':begin
        fix_ind = long(strmid(uvalue,3,strlen(uvalue)-3))
        widget_control,info.fixinfo[fix_ind],get_combobox_select=parameters.lumpars[fix_ind].fixed
     end
     ELSE: print,"Unassigned: "+uvalue
  ENDCASE

END

