pro chain

!p.charsize=1.5

res = mrdfits('output.fits',4,/silent)

p1 = res.p
p2 = res.q
n = n_elements(res)              

set_plot,'x'
plot,p1,p2,psym=-2,xstyle=1,ystyle=1,xrange=[0,-7],yrange=[0,6],xtitle='P',ytitle='Q'
xyouts,p1[0]+0.1,p2[0]+0.1,"Start"
xyouts,p1[n-1]+0.1,p2[n-1]+0.1,"End"

end
