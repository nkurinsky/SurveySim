;Makes model templates, modified blackbody B(v)v^beta
;creates a fits file as output
;Current as of 6/22

pro templates

  !p.thick=5
  !x.thick=5
  !y.thick=5
  !p.charthick=4
  !p.charsize=1.5
  !p.multi=[0,3,4]

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
     t_it = 0

     for i=temp_low,temp_high,temp_step do begin
        res = mod_bb(i,freq,beta)
        for copy=0,n_elements(res)-1 do begin
           templates(b_it,t_it,copy)=res[copy]
        endfor
        t_it++
     endfor
     b_it++
  endfor
  
  print,templates(0,0,0),templates(0,0,1),templates(0,0,2)
  gpts = where(templates eq templates(0,0,1))
  print,gpts

  set_plot,'ps'
  !p.multi=[0,3,4]
  device,filename='plots/template_colors.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated
  
  f250 = dblarr(bnum*tnum)
  f350 = dblarr(bnum*tnum)
  f500 = dblarr(bnum*tnum)
  colors = dblarr(bnum*tnum)

  inds = indgen(tnum)
  for i=0,bnum do begin
     colors[i*bnum+inds] = (i*20+20)
  endfor

  print,n_elements(f250)

  for z=0.d0,10d0,1.0d0 do begin
     colors = [0]
     b250 = 250d-6/(1+z)
     b350 = 350d-6/(1+z)
     b500 = 500d-6/(1+z)
     for b=0,bnum-1 do begin
        for t=0,tnum-1 do begin
           f250[b*10+t] = interpol(reform(templates(b,t,*)),true_wave,b250)
           f350[b*10+t] = interpol(reform(templates(b,t,*)),true_wave,b350)
           f500[b*10+t] = interpol(reform(templates(b,t,*)),true_wave,b500)
           colors = [colors,2.0*(t*temp_step+temp_low)]
        endfor
     endfor

     c1 = acolor(f250,f350,b250,b350)
     c2 = acolor(f350,f500,b350,b500)
     colors = colors[where(colors ne 0)]

     if((z eq 0) or (z eq 10)) then begin
        temp = indgen(bnum*tnum)
        temp = temp[where(temp mod 10 eq 9)]
        print,c1[temp]
        print,c2[temp]
     endif

     loadct,0,/silent
     multiplot
     switch z of
        0.d0:
        3.d0:
        6.d0: begin
           plot,/nodata,[-5,3],[-5,3],xstyle=1,ystyle=1,ytitle=textoidl('\alpha_{350}^{500}')
           break
        end
        1.d0:
        2.d0:
        4.d0:
        5.d0: 
        7.d0: begin
           plot,/nodata,[-5,3],[-5,3],xstyle=1,ystyle=1
           break
        end
        9.d0: begin
           plot,/nodata,[-5,3],[-5,3],xstyle=1,ystyle=1,xtitle=textoidl('\alpha_{250}^{350}'),ytitle=textoidl('\alpha_{350}^{500}')
           break
        end
        8.d0:
        10.d0: begin
           plot,/nodata,[-5,3],[-5,3],xstyle=1,ystyle=1,xtitle=textoidl('\alpha_{250}^{350}')
           break
        end
        else: print,'not found: ',z
     endswitch
     
     loadct,39,/silent
     for i=0,n_elements(c1)-1 do begin
        oplot,[c1[i],c1[i]],[c2[i],c2[i]],color=colors[i],psym=2,symsize=0.25
     endfor
     xyouts,0,2,textoidl('z = '+strmid(string(long(z)),10,2))
  endfor

  multiplot
  plot,/nodata,xstyle=4,ystyle=4,[-5,3],[-5,3]
  legend,['10 K','<->','100 K'],textcolors=[20,0,200],/horizontal,pos=[-4.5,-1],spacing=0,charsize=1.25

  device,/close

  wstruct = {wavelength:true_wave,models:templates}
  elems = n_elements(true_wave)
  wave_sep = (true_wave[elems-1]-true_wave[0])/(elems-1)

  ;; sxaddpar, hdr, 'DATE', systime(),'Date of creation'
  ;; sxaddpar, hdr, 'MODEL','B(T)v','Template Model Used'
  ;; sxaddpar, hdr, 'WAVE_MIN',true_wave[0],'Domain Lower Bound'
  ;; sxaddpar, hdr, 'WAVE_MAX',true_wave[elems-1],'Domain Upper Bound'
  ;; sxaddpar, hdr, 'WAVE_SEP',wave_sep,'Domain Step Size'
  ;; sxaddpar, hdr, 'BETA_MAX',beta_high,'Beta Model Upper Limit'
  ;; sxaddpar, hdr, 'BETA_MIN',beta_low,'Beta Model Lower Limit'
  ;; sxaddpar, hdr, 'BETA_SEP',beta_step,'Beta Model Step Size'
  ;; sxaddpar, hdr, 'TEMP_MAX',temp_high,'Temperature Model Upper Limit'
  ;; sxaddpar, hdr, 'TEMP_MIN',temp_low,'Temperature Model Lower Limit'
  ;; sxaddpar, hdr, 'TEMP_SEP',temp_step,'Temperature Model Step Size'
  ;; mwrfits,templates,'model.fits',hdr,/create
  ;; sxaddpar, hdr2, 'NOTE',"Same as image"
  ;; mwrfits,wstruct,'model.fits',hdr2

  wave = true_wave
  save,filename='templates.save',templates,wave

  ;;write_fits_sim,templates,true_wave,2,[beta_step,temp_step],[beta_low,temp_low],[beta_high,temp_high],model_form='B_v(T)'

end
