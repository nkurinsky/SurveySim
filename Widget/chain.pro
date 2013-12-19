pro chain

!p.charsize=1.5

device,decomposed=0
WINDOW, 0, TITLE='Parameter Space'
res0 = mrdfits('output.fits',0,head,/silent)
print,head
res = mrdfits('output.fits',4,/silent)

p = [res.p00,res.p01,res.p02,res.p03,res.p04]
q = [res.q00,res.q01,res.q02,res.q03,res.q04]
chis = [res.chisq0,res.chisq1,res.chisq2,res.chisq3,res.chisq4]

prange=[min(p),max(p)]
qrange=[min(q),max(q)]
dp=0.1
dq=0.1

set_plot,'x'
plot,prange,qrange,xstyle=1,ystyle=1,/nodata,xtitle='P',ytitle='Q'
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

window,1,title="Chi Square for Chains"
loadct,1
plot,[0,n_elements(p1)],[100,100000],/ylog,xstyle=1,ystyle=1,/nodata,xtitle='Run Number',ytitle='Chisq'
loadct,39
oplot,res.chisq0,color=40
oplot,res.chisq1,color=80
oplot,res.chisq2,color=120
oplot,res.chisq3,color=160
oplot,res.chisq4,color=200

window,2,title="Chi Square for P"
loadct,1
plot,prange,[100,100000],/ylog,xstyle=1,ystyle=1,/nodata,xtitle='P',ytitle='Chisq'
loadct,39

for i=double(prange[0]),double(prange[1]),dp do begin
   print,i
   gpts = where((p gt i) and (p lt i+dp))
   if(gpts[0] ne -1) then begin
      mtemp = mean(chis[gpts])
      mlow = min(chis[gpts])
      mhigh = max(chis[gpts])
      oplot,[i+dp/2],[mtemp],psym=2
      oploterror,[i+dp/2],[mtemp],[mtemp-mlow],/lobar
      oploterror,[i+dp/2],[mtemp],[mhigh-mtemp],/hibar
   endif
endfor  

oplot,res.p00,res.chisq0,color=40,psym=3
oplot,res.p01,res.chisq1,color=80,psym=3
oplot,res.p02,res.chisq2,color=120,psym=3
oplot,res.p03,res.chisq3,color=160,psym=3
oplot,res.p04,res.chisq4,color=200,psym=3

window,3,title="Chi Square for Q"
loadct,1
plot,qrange,[100,100000],/ylog,xstyle=1,ystyle=1,/nodata,xtitle='Q',ytitle='Chisq'
loadct,39
for i=double(qrange[0]),double(qrange[1]),dq do begin
   print,i
   gpts = where((q gt i) and (q lt i+dq))
   if(gpts[0] ne -1) then begin
      mtemp = mean(chis[gpts])
      mlow = min(chis[gpts])
      mhigh = max(chis[gpts])
      oplot,[i+dq/2],[mtemp],psym=2
      oploterror,[i+dq/2],[mtemp],[mtemp-mlow],/lobar
      oploterror,[i+dq/2],[mtemp],[mhigh-mtemp],/hibar
   endif
endfor   

oplot,res.q00,res.chisq0,color=40,psym=3
oplot,res.q01,res.chisq1,color=80,psym=3
oplot,res.q02,res.chisq2,color=120,psym=3
oplot,res.q03,res.chisq3,color=160,psym=3
oplot,res.q04,res.chisq4,color=200,psym=3

end
