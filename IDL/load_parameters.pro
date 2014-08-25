function load_parameters,saveFile=saveFile
  
  if( not keyword_set(saveFile)) then begin
     spawn,'echo $HOME',home
     saveFile=home+'/.SurveySimParams.save'
  endif
  
  if(file_test(saveFile)) then begin
     restore,saveFile
  endif else begin
     
     filter_properties = {FilterProperties,fmin:25.0,ferr:6.2}
     filter = {Filter,filter_id:0,properties:filter_properties}
     filters = replicate(filter,3)
     filters[1].properties = {fmin:20.0,ferr:5.8}
     filters[2].properties = {fmin:15.0,ferr:6.2}
     
                                ;luminosity function parameter initialization
     lpars = {ParameterProperties,value:-2.2,min:-2.4,max:-2.0}
     lumstruct = {Parameter,pars:lpars,fixed:1,name:"PHI0"}
     lumpars = replicate(lumstruct,8)
     
     names = ["PHI0","L0","ALPHA","BETA","P","Q","ZCUT","CEXP"]
     values = [-2.2,10.14,0.5,3.0,-4.5,4.5,2.0,0.0]
     mins = [-2.4,10.0,0.0,0.0,-8.0,-2.0,0.0,0.0]
     maxes = [-2.0,10.3,1.0,5.0,-2.0,8.0,4.0,4.0]
     
     for i=0,n_elements(names)-1 do begin
        lumpars[i].name = names[i]
        lumpars[i].pars.value = values[i]
        lumpars[i].pars.min = mins[i]
        lumpars[i].pars.max = maxes[i]
     endfor
     
                                ;default p and q to fitted parameters
     lumpars[4].fixed = 0
     lumpars[5].fixed = 0
     
     surveyData = {SurveyProperties,area:10.0,zmin:0.0001,zmax:4.0,dz:0.1,runs:1.e4}
     
     CD, Current=thisdir
     files = {FileNames, $
              savefile: saveFile,$
              ofile:'/usr/local/surveysim/obs/spire_fls_dr2.fits', $
              mfile:thisdir+'/model.fits', $
              sedfile:'/usr/local/surveysim/templates/sf_templates.fits', $
              oname:thisdir+'/output.fits' }
     
     msettings = {MCMCSettings,$
                  nchain:5,$
                  tmax:20.0,$
                  acceptpct:0.25, $
                  pct_range:0.05, $
                  conv_conf:0.05, $
                  conv_rmax:1.05, $
                  conv_step:20, $
                  burn_step:10, $
                  burn_ratio:10}
     
     parameters = {SurveySimParameters, $
                   filters:filters,$
                   lumpars:lumpars,$
                   surveyData:surveyData,$
                   files:files,$
                   msettings:msettings,$
                   print:1}
  endelse

  return,parameters
  
end
  
