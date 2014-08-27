function binsize,data,min,max
  stats = moment(data,maxmoment=2)
  datanum = n_elements(data)
  range = max-min
  
  db = 3.49*sqrt(stats[1])/(datanum^0.333333)
  nbins = floor(range/db)
  db = range/float(nbins)
  return,db
end
