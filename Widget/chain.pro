pro chain

!p.charsize=1.5

device,decomposed=0
res = mrdfits('output.fits',4,/silent)

p1 = res.p
p2 = res.q
n = n_elements(res)              

set_plot,'x'
plot,[0,-7],[0,6],xstyle=1,ystyle=1,/nodata,xtitle='P',ytitle='Q'

loadct,39

oplot,p1,p2,psym=2
xyouts,p1[0]+0.1,p2[0]+0.1,"Start"
xyouts,p1[n-1]+0.1,p2[n-1]+0.1,"End"

end
