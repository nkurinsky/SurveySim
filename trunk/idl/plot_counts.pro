pro plot_counts,filename,bandnum

  dists = mrdfits(filename,3,head,/silent)
  count_dists = mrdfits(filename,6,head,/silent)

  alpha = 0.159                 ;one std deviation
  plusfrac = (1.0-alpha)
  minusfrac = alpha

  Case bandnum of
     1: begin
        x = dists.s1
        ybest = dists.mod_dnds1
        yobs = dists.obs_dnds1
     end
     2: begin
        x = dists.s2
        ybest = dists.mod_dnds2
        yobs = dists.obs_dnds2
     end
     3: begin
        x = dists.s3
        ybest = dists.mod_dnds3
        yobs = dists.obs_dnds3
     end
     else: print,"Invalid band number "+strtrim(string(bandnum),1)
  endcase

  if((bandnum gt 0) and (bandnum lt 4)) then begin
     c = count_dists.(bandnum-1)
     csize = n_elements(c[*,0])
     cmedian = make_array(csize,value=0.0)
     cplus = make_array(csize,value=0.0)
     cminus = make_array(csize,value=0.0)

     for i=0,csize-1 do begin
        dnds = c[csize-i-1,*]
        dnds = dnds[sort(dnds)]
        dnds = dnds[where(dnds gt 0)]
        pi = plusfrac*n_elements(dnds)
        mi = minusfrac*n_elements(dnds)
        cmedian[i] = median([dnds])
        cplus[i] = dnds[pi]
        cminus[i] = dnds[mi]
     endfor

     bestpts = where(ybest gt 0)
     xbest = x[bestpts]/1.d3
     ybest = ybest[bestpts]

     obspts = where(yobs gt 0)
     xobs = x[obspts]/1.d3
     yobs = yobs[obspts]

     minpts = where(cminus gt 0)
     xstats = x[minpts]/1.d3
     yminus = cminus[minpts]
     yplus = cplus[minpts]
     ymedian = cmedian[minpts]
     
     plot,xobs,yobs, $
          /xlog,/ylog, $
          psym=2,symsize=2, $
          xtitle=TeXtoIDL('F_{'+strtrim(string(bandnum),1)+'}[Jy]'), $
          ytitle=TeXtoIDL('(dN/dS)S^{2.5} [gal ster^{-1} J^{1.5}]'), $
          title='Band '+strtrim(string(bandnum),1)+' Counts'

     oplot,xbest,ybest,psym=4,symsize=2
     oplot,xstats,ymedian,linestyle=2
     oplot,xstats,yplus,linestyle=1
     oplot,xstats,yminus,linestyle=1
  endif

end
