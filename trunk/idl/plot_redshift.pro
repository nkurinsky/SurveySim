pro plot_redshift,filename,zmin,zmax,dz
  
  dists = mrdfits(filename,3,head,/silent)
  gpts = where(dists.f3 gt 0)
  z = dists[gpts].z
  
  h = histogram(z,binsize=dz,locations=xh,min=zmin,max=zmax)
  plot,xh,h,psym=10,xrange=[zmin,zmax],xstyle=1,xtitle='z',ytitle='dN/dz',title='Redshift Distribution'

end
