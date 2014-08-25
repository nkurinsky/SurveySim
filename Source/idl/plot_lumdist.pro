pro plot_lumdist,filename
  
  dists = mrdfits(filename,3,head,/silent)
  gpts = where(dists.f3 gt 0)
  lum = dists[gpts].lum
  f = alog10(lum)

  h = histogram(f,locations=xh)
  plot,xh,h,psym=10,/xlog,xstyle=1,ystyle=0,title="Luminosity Distribution",xtitle="L",ytitle="dN/dL"
  
end
