;Makes model templates, modified blackbody B(v)v^beta
;creates a fits file as output
;Current as of 6/19/2010

pro templates

  !p.thick=5
  !x.thick=5
  !y.thick=5
  !p.charthick=4
  !p.charsize=0.75
  !p.multi=[0,3,4]

  set_plot,'ps'
  device,filename='plots/template.eps',xsize=14,ysize=12,/inches,/times,/color,/encapsulated
  loadct,39,/silent
  wave = (findgen(991)+10)*1.d4
  true_wave = wave*1d-10
  freq = 2.9979d8/true_wave
  plot_wave = wave*1d-4

  beta_high = 2.d0
  beta_low = 1.d0
  beta_step = 0.1d0
  temp_high = 100.d0
  temp_low = 10.d0
  temp_step = 10.d0

  bnum = ((beta_high+beta_step-beta_low)/beta_step)
  tnum = ((temp_high-temp_low+temp_step)/temp_step)
  templates = dblarr(bnum,tnum,n_elements(wave))

  b_it = 0
  for beta=beta_low,beta_high+beta_step,beta_step do begin
     print,beta
     lowest = mod_bb(10,2.9979e12,beta)
     highest = max(mod_bb(100,freq,beta))
     plot,/nodata,xstyle=1,ystyle=1,[1000,10],[lowest,10^(alog10(lowest)+8)],/xlog,/ylog,xtitle=textoidl('\lambda [\mum]'),ytitle=textoidl('B_{\nu}(T[K])\nu^{'+strmid(strcompress(string(beta+0.005),/remove_all),0,4)+'}'),charsize=1.5
     t_it = 0

     for i=temp_low,temp_high,temp_step do begin
        ;res = blackbody(i,freq,/freq)*freq^beta
        res = mod_bb(i,freq,beta)
        oplot,plot_wave,res,color=2.0*i
        if((i eq temp_low) and (beta eq beta_low)) then print,res[0],res[100],res[900]
        for copy=0,n_elements(res)-1 do begin
           templates(b_it,t_it,copy)=res[copy]
        endfor
        t_it++
     endfor
     b_it++
  endfor
  
  print,templates(0,0,0),templates(0,0,1),templates(0,0,2)
  print,templates(1,0,0),templates(1,0,1),templates(1,0,2)

  plot,/nodata,xstyle=4,ystyle=4,[0,1],[0,1]
  legend,['Temp:',' 10 K',' 20 K',' 30 K',' 40 K',' 50 K',' 60 K',' 70 K',' 80 K',' 90 K','100 K'],textcolors=findgen(11)*20,charsize=1.25

  device,/close

  !p.multi=[0,3,4]
  device,filename='plots/template_colors.eps',xsize=14,ysize=12,/inches,/times,/color,/encapsulated
  
  f250 = dblarr(bnum*tnum)
  f350 = dblarr(bnum*tnum)
  f500 = dblarr(bnum*tnum)
  print,n_elements(f250)

  for z=0.d0,10d0,1.0d0 do begin
     b250 = replicate(250d-6,n_elements(f250))/(1+z)
     b350 = replicate(350d-6,n_elements(f250))/(1+z)
     b500 = replicate(500d-6,n_elements(f250))/(1+z)
     for b=0,bnum-1 do begin
        for t=0,tnum-1 do begin
           f250[b*10+t] = interpol(reform(templates(b,t,*)),true_wave,b250[0])
           f350[b*10+t] = interpol(reform(templates(b,t,*)),true_wave,b350[0])
           f500[b*10+t] = interpol(reform(templates(b,t,*)),true_wave,b500[0])
        endfor
     endfor

     c1 = acolor(f250,f350,b250,b350)
     c2 = acolor(f350,f500,b350,b500)
     
     plot,c1,c2,psym=2,xrange=[-5,3],yrange=[-5,3],xstyle=1,ystyle=1,charsize=1.5,xtitle=textoidl('\alpha_{250}^{350}'),ytitle=textoidl('\alpha_{350}^{500}')
     xyouts,1,2,textoidl('z = '+strmid(strcompress(string(z),/remove_all),0,3)),charsize=1.5
  endfor

  device,/close

  wstruct = {wavelength:true_wave,models:templates}
  elems = n_elements(true_wave)
  wave_sep = (true_wave[elems-1]-true_wave[0])/(elems-1)

  sxaddpar, hdr, 'DATE', systime(),'Date of creation'
  sxaddpar, hdr, 'MODEL','B(T)v','Template Model Used'
  sxaddpar, hdr, 'WAVE_MIN',true_wave[0],'Domain Lower Bound'
  sxaddpar, hdr, 'WAVE_MAX',true_wave[elems-1],'Domain Upper Bound'
  sxaddpar, hdr, 'WAVE_SEP',wave_sep,'Domain Step Size'
  sxaddpar, hdr, 'BETA_MAX',beta_high,'Beta Model Upper Limit'
  sxaddpar, hdr, 'BETA_MIN',beta_low,'Beta Model Lower Limit'
  sxaddpar, hdr, 'BETA_SEP',beta_step,'Beta Model Step Size'
  sxaddpar, hdr, 'TEMP_MAX',temp_high,'Temperature Model Upper Limit'
  sxaddpar, hdr, 'TEMP_MIN',temp_low,'Temperature Model Lower Limit'
  sxaddpar, hdr, 'TEMP_SEP',temp_step,'Temperature Model Step Size'
  mwrfits,templates,'model.fits',hdr,/create
  sxaddpar, hdr2, 'NOTE',"Same as image"
  mwrfits,wstruct,'model.fits',hdr2

end
