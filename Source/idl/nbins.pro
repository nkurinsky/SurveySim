function nbins,data,min,max
  stats = moment(data,maxmoment=2)
  datanum = n_elements(data)
  db = 3.49*sqrt(stats[1])/(datanum^0.333333)
  nbins = floor((max-min)/db)
  return,long(nbins)
end
