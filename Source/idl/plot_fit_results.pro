pro plot_fit_results,filename
  
  res = mrdfits(filename,4,/silent)
  alltags = tag_names(res)
  tags = alltags[where(strmatch(alltags,'*0',/FOLD_CASE) eq 1)]
  
  for i=0,n_elements(tags)-1 do begin
     tags[i] = strmid(tags[i],0,strlen(tags[i])-1)
  endfor
  tags = tags[where((tags ne 'CHISQ') and (tags ne 'ACPT'))]
  print,'Fitted Variables: ',tags

  dim = n_elements(tags)
  nbins = 20.0
  
  cpts = where(strmatch(alltags,'CHISQ*',/FOLD_CASE) eq 1)
  for ci=0,n_elements(cpts)-1 do begin
     if(ci eq 0) then begin
        chis = res.(cpts[ci])
     endif else begin
        chis = [chis,res.(cpts[ci])]
     endelse
  endfor
  
  apts = where(strmatch(alltags,'ACPT*',/FOLD_CASE) eq 1)
  accept = make_array(n_elements(apts),n_elements(res.(cpts[0])),value=0.0)
  for ai=0,n_elements(apts)-1 do begin
     if(ai eq 0) then begin
        accept = res.(apts[ai])
     endif else begin
        accept = [accept,res.(apts[ai])]
     endelse
  endfor

  chi_med = median(chis)
  print,"median: ",chi_med
  chistpts = where(chis lt chi_med)

  chi_min = alog(min(chis))
  chi_max = alog(max(chis))
  color_scale = 256/n_elements(chistpts)

  multiplot,[dim,dim],/init,/rowmajor,mTitle="MCMC Fitting Results",gap=0.005
  cgText, 0.6, 0.9, Alignment=0, /Normal, 'Fitting Results:', Charsize=1.25

  for x=0,dim-1 do begin
     
     ppts = where(strmatch(alltags,tags[x]+'*',/FOLD_CASE) eq 1)
     for pi=0,n_elements(ppts)-1 do begin
        if(pi eq 0) then begin
           p = res.(ppts[pi])
        endif else begin
           p = [p,res.(ppts[pi])]
        endelse
     endfor
     prange = [min(p),max(p)]

     for y=0,dim-1 do begin
        multiplot

        if (x eq y) then begin  ;make histogram

           pgood = p[chistpts]
           stats = moment(pgood,maxmoment=2)
           print,"Variable: ",tags[x],", Mean: ",stats[0]," Variance: ",stats[1]
           cgText, 0.6, 0.87-0.03*x, Alignment=0, /Normal, Charsize=1.25,textoidl(tags[x]+": "+strcompress(string(stats[0],format='(D0.3)'))+"\pm"+strcompress(string(stats[1],format='(D0.3)')))
           h = histogram(pgood,locations=xh, $
                         binsize=3.49*sqrt(stats[1])/(n_elements(pgood)^0.333333))
           
           h = float(h)/float(total(h))
           loadct,1,/silent
           
           xt = ''
           yt = ''
           if(x eq 0) then yt=tags[y]
           if(y eq dim-1) then xt=tags[x]
           plot,xh,h,psym=10,xrange=prange,yrange=[0,1],xtitle=xt,ytitle=yt,xstyle=1,ystyle=1

        endif else if (x lt y) then begin ;make liklihood space

           qpts = where(strmatch(alltags,tags[y]+'*',/FOLD_CASE) eq 1)
           for qi=0,n_elements(qpts)-1 do begin
              if(qi eq 0) then begin
                 q = res.(qpts[qi])
              endif else begin
                 q = [q,res.(qpts[qi])]
              endelse
           endfor
           qrange = [min(q),max(q)]
           qgood = q[chistpts]
           
           dp=3.49*sqrt(stats[1])/(n_elements(pgood)^0.333333)
           dq=3.49*sqrt(stats[1])/(n_elements(qgood)^0.333333)

           xt = ''
           yt = ''
           if(x eq 0) then yt=tags[y]
           if(y eq dim-1) then xt=tags[x]
           
           loadct,1,/silent
           plot,prange,qrange,/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=xt,ytitle=yt
           loadct,39,/silent

           for i=prange[0],prange[1]-dp/2,dp do begin
              for j=qrange[0],qrange[1]-dq/2,dq do begin

                 gpts = where((pgood gt i) and (pgood le i+dp) and (qgood gt j) and (qgood lt j+dq))
                 
                 if(gpts[0] ne -1) then begin
                    xfill = [i,i,i+dp,i+dp]
                    yfill = [j,j+dq,j+dq,j]
                    polyfill,xfill,yfill,color=color_scale*n_elements(gpts)
                 endif
              endfor
           endfor

        endif 
     endfor
  endfor
  
  multiplot,/reset

end
