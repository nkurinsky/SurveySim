function findex,u,v
   nu=n_elements(u)
   nv=n_elements(v)

   us=u-shift(u,+1)
   us=us(1:*)
   umx=max(us,min=umn)
   if umx gt 0 and umn lt 0 then message,'u must be monotonic'
   if umx gt 0 then inc=1 else inc=0

   maxcomp=fix(alog(float(nu))/alog(2.)+.5) 

   ; maxcomp = maximum number of binary search iteratios

   jlim=lonarr(2,nv)
   jlim(0,*)=0          ; array of lower limits
   jlim(1,*)=nu-1       ; array of upper limits

   iter=0
   repeat begin
     jj=(jlim(0,*)+jlim(1,*))/2
     ii=where(v ge u(jj),n) & if n gt 0 then jlim(1-inc,ii)=jj(ii)
     ii=where(v lt u(jj),n) & if n gt 0 then jlim(inc,ii)=jj(ii)
     jdif=max(jlim(1,*)-jlim(0,*))
     if iter gt maxcomp then begin
       print,maxcomp,iter, jdif
       message,'binary search failed'
     endif
     iter=iter+1
   endrep until jdif eq 1 

   w=v-v
   w(*)=(v-u(jlim(0,*)))/(u(jlim(0,*)+1)-u(jlim(0,*)))+jlim(0,*)

return,w
end
