function get_parameters,saveFile
  
  if(file_test(saveFile)) then begin
     restore,saveFile
  endif else begin
     
  filter_properties = {fmin:25.0,ferr:6.2}
  filter = {filter_id:1,properties:filter_properties}
  filters = replicate(filter,3)
  
  filters[1].filter_id=2
  filters[1].properties = {fmin:20.0,ferr:5.8}
  
  filters[2].filter_id=3
  filters[2].properties = {fmin:15.0,ferr:6.2}

  ;luminosity function parameter initialization
  lpars = {value:-2.2,min:-2.4,max:-2.0}
  lumstruct = {pars:lpars,fixed:1,name="PHI0"}
  lumpars = replicate(lumstruct,8)

  names = ["PHI0","L0","ALPHA","BETA","P","Q","ZCUT","CEXP"]
  values = [-2.2,10.14,0.5,3.0,-4.5,4.5,2.0,cexp]
  mins = [-2.4,10.0,0.0,0.0,-8.0,-2.0,0.0,0.0]
  maxes = [-2.0,10.3,1.0,5.0,-2.0,8.0,4.0,4.0]
  
  for i,n_element(names-1) do begin
     lumpars[i].name = names[i]
     lumpars[i].pars.value = values[i]
     lumpars[i].pars.min = mins[i]
     lumpars[i].pars.max = maxes[i]
  endfor

  ;default p and q to fitted parameters
  lumpars[4].fix = 0
  lumpars[5].fix = 0
  
  surveyData = {area:10.0,zmin:0.0001,zmax:4.0,dz:0.1,runs:1.e4}
  
  CD, Current=thisdir
  files = {ofile:'/usr/local/surveysim/obs/spire_fls_dr2.fits', $
           mfile:thisdir+'/model.fits', $
           sedfile:'/usr/local/surveysim/templates/sf_templates.fits', $
           oname:thisdir+'/output.fits' }
  
  msettings = {$
              nchain:5,$
              tmax:20.0,$
              acceptpct:0.25, $
              pct_range:0.05, $
              conv_conf:0.05, $
              conv_rmax:1.05, $
              conv_step:20, $
              burn_step:10, $
              burn_ratio:10}

  parameters = {filters:filters,$
                lumpars:lumpars,$
                surveyData:surveyData,$
                files:files,$
                msettings:msettings,$
                print:1}
  
  return,parameters
end
