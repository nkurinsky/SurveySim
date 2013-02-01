;sets defaults and runs the functions that make histograms of first
;model, then observed, and finally the comparison
;Current as of 6/22

pro run_comp

  !p.multi=[0,1,0]
  
  min = -5.0
  max = 10.0
  bin = 0.5

  model_color,min,max,bin
  obs_color,min,max,bin
  compare,min,max,bin

end
