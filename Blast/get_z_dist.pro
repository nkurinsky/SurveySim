function get_z_dist

  near = randomn(seed,2e3)*4.0
  near = abs(near)
  mid = randomn(seed,6e3)
  mid -= min(mid)
  mid*=2.5/max(mid)
  redshift = randomn(seed,1e4)
  redshift -= min(redshift)
  redshift *= 6.0/max(redshift)
  redshift = [redshift,near,mid]

  gpts = where((redshift le 9.0) and (redshift ge 0))
  redshift = redshift[gpts]

  print,'Redshift:'
  print,'Mean:',mean(redshift)
  print,'Low Percentage:',float(n_elements(where(redshift lt 1.2)))/float(n_elements(redshift))

  return,redshift

end
