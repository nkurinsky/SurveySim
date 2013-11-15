pro chain

!p.charsize=1.5

device,decomposed=0
res = mrdfits('output.fits',4,/silent)

set_plot,'x'
plot,[-6,2],[0,8],xstyle=1,ystyle=1,/nodata,xtitle='P',ytitle='Q'
n = n_elements(res) 
loadct,39

p1 = res.p00
p2 = res.q00
oplot,p1,p2,psym=3,color=40
xyouts,p1[0]+0.1,p2[0]+0.1,"Start"
xyouts,p1[n-1]+0.1,p2[n-1]+0.1,"End"

p1 = res.p01
p2 = res.q01
oplot,p1,p2,psym=3,color=80
xyouts,p1[0]+0.1,p2[0]+0.1,"Start"
xyouts,p1[n-1]+0.1,p2[n-1]+0.1,"End"

p1 = res.p02
p2 = res.q02
oplot,p1,p2,psym=3,color=120
xyouts,p1[0]+0.1,p2[0]+0.1,"Start"
xyouts,p1[n-1]+0.1,p2[n-1]+0.1,"End"

p1 = res.p03
p2 = res.q03
oplot,p1,p2,psym=3,color=160
xyouts,p1[0]+0.1,p2[0]+0.1,"Start"
xyouts,p1[n-1]+0.1,p2[n-1]+0.1,"End"

p1 = res.p04
p2 = res.q04
oplot,p1,p2,psym=3,color=200
xyouts,p1[0]+0.1,p2[0]+0.1,"Start"
xyouts,p1[n-1]+0.1,p2[n-1]+0.1,"End"

end
