pro plot_counts,filename,bandnum

  dists = mrdfits(filename,3,head,/silent)
  counts = compute_counts(filename)
  
  Case bandnum of
     1: begin
        xd1 = dists.s1
        yd1 = dists.mod_dnds1
        gpts = where(yd1 gt 0)
        xd1 = xd1[gpts]/1.d3
        yd1 = yd1[gpts]
        
        if( file_test('counts_clements10.dat')) then begin
           readcol,'counts_clements10.dat',skipline=2,numline=16,flux,nbin,corr,int_counts,int_err,diff_counts,diff_err,/silent
           flux /= 1.d3
           xrange=[min([xd1,flux]),max([xd1,flux])]
           gpts = where(counts.band1.minus gt 0)
           yrange=[min([yd1,diff_counts,counts.band1.minus[gpts]])/1.2,max([yd1,diff_counts,counts.band1.plus])*1.2]
           plot,flux,diff_counts,psym=1,symsize=2,xtitle=TeXtoIDL('F_{250}[Jy]'),ytitle=TeXtoIDL('(dN/dS)S^{2.5} [gal ster^{-1} J^{1.5}]'),/xlog,/ylog,title='Band 1 Counts',yrange=yrange,xrange=xrange,ystyle=1,xstyle=1
           oploterr,flux,diff_counts,diff_err
           oplot,xd1,yd1,psym=2
        endif else begin
           print,'Error: File "counts_clements10.dat" not found'
           plot,xd1,yd1,psym=2,symsize=2,xtitle=TeXtoIDL('F_{250}[Jy]'),ytitle=TeXtoIDL('(dN/dS)S^{2.5} [gal ster^{-1} J^{1.5}]'),/xlog,/ylog,title='Band 1 Counts'
        endelse
        oplot,xd1,counts.band1.mean,linestyle=2
        oplot,xd1,counts.band1.plus,linestyle=1
        oplot,xd1,counts.band1.minus,linestyle=1
     end
     2: begin
        xd2 = dists.s2
        yd2 = dists.mod_dnds2
        gpts = where(yd2 gt 0)
        xd2 = xd2[gpts]/1.d3
        yd2 = yd2[gpts]
        
        if( file_test('counts_clements10.dat')) then begin
           readcol,'counts_clements10.dat',skipline=19,numline=13,flux,nbin,corr,int_counts,int_err,diff_counts,diff_err,/silent
           flux /= 1.d3
           xrange=[min([xd2,flux]),max([xd2,flux])]
           gpts= where(counts.band2.minus gt 0)
           yrange=[min([yd2,diff_counts,counts.band2.minus[gpts]])/1.2,max([yd2,diff_counts,counts.band2.plus])*1.2]
           plot,flux,diff_counts,psym=1,symsize=2,xtitle=TeXtoIDL('F_{250}[Jy]'),ytitle=TeXtoIDL('(dN/dS)S^{2.5} [gal ster^{-1} J^{1.5}]'),/xlog,/ylog,title='Band 2 Counts',yrange=yrange,xrange=xrange,ystyle=1,xstyle=1
           oploterr,flux,diff_counts,diff_err
           oplot,xd2,yd2,psym=2
        endif else begin
           plot,xd2,yd2,psym=2,symsize=2,xtitle=TeXtoIDL('F_{250}[Jy]'),ytitle=TeXtoIDL('(dN/dS)S^{2.5} [gal ster^{-1} J^{1.5}]'),/xlog,/ylog,title='Band 2 Counts'
        endelse
        oplot,xd2,counts.band2.mean,linestyle=2
        oplot,xd2,counts.band2.plus,linestyle=1
        oplot,xd2,counts.band2.minus,linestyle=1
        
     end
     3: begin
        xd3 = dists.s3
        yd3 = dists.mod_dnds3
        gpts = where(yd3 gt 0)
        xd3 = xd3[gpts]/1.d3
        yd3 = yd3[gpts]
        
        if( file_test('counts_clements10.dat')) then begin
           readcol,'counts_clements10.dat',skipline=33,numline=10,flux,nbin,corr,int_counts,int_err,diff_counts,diff_err,/silent
           flux /= 1.d3
           xrange=[min([xd3,flux]),max([xd3,flux])]
           gpts= where(counts.band3.minus gt 0)
           yrange=[min([yd3,diff_counts,counts.band3.minus[gpts]])/1.2,max([yd3,diff_counts,counts.band3.plus])*1.2]
           plot,flux,diff_counts,psym=1,symsize=2,xtitle=TeXtoIDL('F_{250}[Jy]'),ytitle=TeXtoIDL('(dN/dS)S^{2.5} [gal ster^{-1} J^{1.5}]'),/xlog,/ylog,title='Band 3 Counts',yrange=yrange,xrange=xrange,ystyle=1,xstyle=1
           oploterr,flux,diff_counts,diff_err
           if(gpts[0] ne -1) then oplot,xd3,yd3,psym=2
        endif else if(gpts[0] ne -1) then begin
           plot,xd3,yd3,psym=2,symsize=2,xtitle=TeXtoIDL('F_{250}[Jy]'),ytitle=TeXtoIDL('(dN/dS)S^{2.5} [gal ster^{-1} J^{1.5}]'),/xlog,/ylog,title='Band 3 Counts'
        endif
        
        if(gpts[0] ne -1) then begin
           oplot,xd3,counts.band3.mean,linestyle=2
           oplot,xd3,counts.band3.plus,linestyle=1
           oplot,xd3,counts.band3.minus,linestyle=1
        endif

     end
     else: print,"Invalid band number "+strtrim(string(bandnum),1)
     
  endcase
end
