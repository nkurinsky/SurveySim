;==================================================================
; Writen by Noah Kurinsky, version recent as of 2/15/13
; Edited by Anna Sajina, December 2012
; certain files in the same directory are required for proper
; function of this widget:
;
;==================================================================

pro simulation

  COMMON simulation_com,info,ldata,ldat0,sdat,cdat,bands

;==================================================================
;the INFO structure holds the key widget control parameters as well as 
;the basic simulation settings
;------------------------------------------------------------------
  info={$ ;widget settings
       base: 0L, $
       base1: 0L, $
       base2: 0L, $
       base_out: 0L, $
       base1_out: 0L, $
       base2_out: 0L, $
       xsize1:150., $
       ysize1:100., $
       xsize2:150., $
       ysize2:100., $
       magnification:1.00, $  
       draw1:0L, $
       draw2:0L, $
       draw3:0L, $
       win_id1:0L, $
       win_id2:0L, $
       win_id3:0L, $
       ncolors:0L, $
;widget bases
       p_main:0L, $
       obs_table:0L, $
       lum_table:0L, $
       sed_table:0L, $
       sim_table:0L, $
       dbase:0L, $
       button_base:0L, $
;input tables and files
       obsname:0L, $
       sfile:0L, $
       ot:0L, $
       t1:0L, $
       t2:0L, $
       t3:0L, $
;basic simulation settings
       ofile:'observation.save', $
       mfile:'model.save', $
       sedfile: 'sf_templates.fits', $
       oname:'output.fits' }

;====================================================================
; Colors
;--------------------------------------------------------------------
;device, truecolor, depth=24 ;decompose=0
  loadct,3,/silent
  tvlct,red,green,blue,/get
  info.ncolors=!d.TABLE_SIZE
  red[info.ncolors-1]=0B
  green[info.ncolors-1]=255B
  blue[info.ncolors-1]=0B
  tvlct,red,green,blue
  
; Screen size
  size_screen = get_screen_size()
  size_screen=size_screen*0.8
  spectrum_xsize=size_screen(0)*2./3. & spectrum_ysize=size_screen(1)/5.
  
  plotSize = 100.
  info.magnification = size_screen[1]/(plotSize+info.ysize1)*0.9
  
;widget base initialization
  info.base = widget_base(title='Model Setup and Initialization',/column,/align_center) ;main base
  info.base2= widget_base(info.base, /column,/align_center)               ;plot base
  info.p_main = widget_base(info.base,/column,/align_center)                   ;parameter base
  info.obs_table = widget_base(info.p_main,/column,/align_center)
  info.lum_table = widget_base(info.p_main,/column,/align_center)
  info.sed_table = widget_base(info.p_main,/column,/align_center)
  info.sim_table = widget_base(info.p_main,/column,/align_center)
  info.dbase = widget_base(info.p_main,/column,/align_center)    ; base for file dialogs
  info.button_base = widget_base(info.p_main,/row,/align_center) ; base for buttons
  
;output file select
  CD, Current=thisDir
  info.obsname = fsc_fileselect(info.dbase,Directory=thisDir,filename=info.ofile,/NoMaxSize,LabelName='Observation Save File') ;,ObjectRef=obsObject) 
  info.sfile = fsc_fileselect(info.dbase,Directory=thisDir,filename=info.sedfile,ObjectRef=sedObject,/NoMaxSize,LabelName='SED Templates file')
  
;=============================================================
;survey parameter initialization
;-------------------------------------------------------------
;don't want to use the existing params.save as the wavelength
;is in meters, whereas we use microns normally
  if(file_test('params.save')) then begin
     restore,'params.save'
  endif else begin
     band1 = {wave:250.d0,fmin:25.0,ferr:6.2}
     band2 = {wave:350.d0,fmin:20.0,ferr:5.8}
     band3 = {wave:500.d0,fmin:15.0,ferr:6.2}
     bands = [band1,band2,band3]

     ;luminosity function parameter initialization
     ldat0={phi0:-2.2,lo:10.14,alpha:0.5,beta:3.0,p:6.7,q:3.5}
     ldat1={phi0:1.0,lo:1.0,alpha:1.0,beta:1.0,p:0.0,q:0.0}
     ldat2={phi0:-5.0,lo:9.0,alpha:0.0,beta:0.0,p:0.,q:0.0}
     ldat3={phi0:5.0,lo:11.0,alpha:2.0,beta:5.0,p:10.0,q:10.0}
     ldat4={phi0:0.0,lo:0.0,alpha:0.0,beta:0.0,p:0.1,q:0.1}
     ldata=[ldat0,ldat1,ldat2,ldat3,ldat4]
     ;the fixed values: =1 if held fixed, =0 if variable

     cdat = {a0:0.0,fixed:0.0,amin:0.0,amax:1.0,sigma:0.1}
     sdat = {area:10.0,zmin:0.0,zmax:5.0,dz:0.1,runs:1.e3}
  endelse

  bname = ["Band 1","Band 2","Band 3"]
  ocols = ["Wavelength (um)","Flux limit (mJy)","Standard Error (mJy)"]
  f = ['(e9.2)','(f7.4)','(f7.4)']
  fmt = [[f],[f],[f]]
  
;The Survey properties table
  lo = widget_label(info.obs_table,value="Survey Properties")
  info.ot = widget_table(info.obs_table,value=bands,column_labels=ocols,row_labels=bname,uvalue='ot',/editable,alignment=1,column_widths=[100,100,150],format=fmt,scr_xsize=425,scr_ysize=95)

  f2a = ['(f5.2)','(f5.2)','(f5.2)','(f5.2)','(f5.2)','(f5.2)']
  f2b = ['(i)','(i)','(i)','(i)','(i)','(i)']
  fmt2 = [[f2a],[f2b],[f2a],[f2a],[f2a]]
  fmt3 = ['(f5.2)','(i)','(f5.2)','(f5.2)','(f5.2)']
  lrows=["Initial","Fixed","Min","Max","Sigma"]

;Luminosity Function Parameters
  l1 = widget_label(info.lum_table,value="Luminosity Function Parameters")
  info.t1 = widget_table(info.lum_table,value=ldata,column_labels=tag_names(ldata),row_labels=lrows,uvalue='t1',/editable,alignment=1,format=fmt2,scr_xsize=472,scr_ysize=132)

  l3 = widget_label(info.sed_table,value="SED Evolution Parameters")
  info.t3 = widget_table(info.sed_table,value=cdat,column_labels=lrows,row_labels=["Color Exp"],uvalue='t3',/editable,alignment=1,format=fmt3,scr_xsize=404,scr_ysize=55)

  tcols = ["Area (sdeg)","Z Min","Z Max","Z Binsize","Run Number"]
;Simulation Parameters
  l2 = widget_label(info.sim_table,value="Simulation Settings")
  info.t2 = widget_table(info.sim_table,value=sdat,column_labels=tcols,/no_row_headers,uvalue='t2',/editable,alignment=1,format=['(f5.2)','(f5.2)','(f5.2)','(f5.2)','(e9.2)'],column_widths=[100,100,100,100,100],scr_xsize=505,scr_ysize=55)
  

;============================================================================
;Show SED templates
;----------------------------------------------------------------------------
info.draw3 = widget_draw(info.base2, xsize=info.magnification*info.xsize2,$
                         ysize=info.magnification*info.ysize1 $
                         ,uvalue="DRAW_WINDOW3",retain=2 $
                         ,/button_events, keyboard_events=1,/tracking_events)

;initialize widget and establish control
widget_control,/realize,info.base,xoffset=0,yoffset=0,Set_UValue=theObject
widget_control, info.draw3, get_value=win_id3 & info.win_id3=win_id3
xmanager,'simulation',info.base,/NO_BLOCK
wset,info.win_id3

;plot SEDs (generalize to any SED template file)
widget_control,info.sfile,get_value=value
templ=mrdfits(info.sedfile,/silent) 

loadct,1,/silent
plot,templ[*,0],templ[*,1],/xlog,/ylog,yrange=[1.d20,1.d28],ystyle=1,xtitle=TeXtoIDL('\lambda [\mum]'),ytitle=TeXtoIDL('L_{\nu} [W/Hz]')
for ipl=1,13 do oplot,templ[*,0],templ[*,ipl+1]

;===========================================================================
;Buttons at bottom 
;---------------------------------------------------------------------------
  run_btn = widget_button(info.button_base,uvalue='go',value='Run Simulation',xsize=100,ysize=25)
  replot_btn = widget_button(info.button_base,uvalue='replot',value='Plot Last Run',xsize=100,ysize=25)
  quit_btn = widget_button(info.button_base,uvalue='quit',value='Quit',xsize=50,ysize=25)
  info_btn = widget_button(info.button_base,uvalue='info',value='Info',xsize=50,ysize=25)

END

;widget event handling routine
PRO simulation_event,ev
  COMMON simulation_com,info,ldata,ldat0,sdat,cdat,bands

  ; get event identifier
  widget_control,ev.id,get_uvalue=uvalue

  CASE uvalue OF
     'go'  : begin ;save settings, intialize FITS file, and pass to C++ fitting routine
        save,ldata,ldat0,bands,sdat,cdat,filename='params.save' ;save parameters
        
        ;make observation FITS file
        widget_control,info.ot,get_value=bvals
        wave = [bvals[0].wave,bvals[1].wave,bvals[2].wave]
        flux_min = [bvals[0].fmin,bvals[1].fmin,bvals[2].fmin]
        flux_err = [bvals[0].ferr,bvals[1].ferr,bvals[2].ferr]

        widget_control,info.obsname,get_value=value
        if (file_test(info.ofile)) then begin
           print,info.ofile
           restore,info.ofile
           values = dblarr(n_elements(wave),n_elements(f1))
        endif else begin
           print,'Error: file "observation.save" not found'
           print,'Simulation attempt unsuccessful'
           break
        endelse

        for i=0,n_elements(f1)-1 do begin
           values(0,i) = f1[i]
           values(1,i) = f2[i]
           values(2,i) = f3[i]
        endfor
        
        sxaddpar, hdr, 'DATE', systime(),'Date of creation'
        sxaddpar, hdr, 'FLUX_NUM',n_elements(f1),'Number of observations included in file'
        sxaddpar, hdr, 'WAVE_1',wave[0],'Wavelength corresponding to first flux'
        sxaddpar, hdr, 'WAVE_2',wave[1],'Wavelength corresponding to second flux'
        sxaddpar, hdr, 'WAVE_3',wave[2],'Wavelength corresponding to third flux'
        sxaddpar, hdr, 'W1_FMIN',flux_min[0],'Wavelength corresponding to first flux'
        sxaddpar, hdr, 'W2_FMIN',flux_min[1],'Wavelength corresponding to second flux'
        sxaddpar, hdr, 'W3_FMIN',flux_min[2],'Wavelength corresponding to third flux'
        sxaddpar, hdr, 'W1_FERR',flux_err[0],'Wavelength corresponding to first flux'
        sxaddpar, hdr, 'W2_FERR',flux_err[1],'Wavelength corresponding to second flux'
        sxaddpar, hdr, 'W3_FERR',flux_err[2],'Wavelength corresponding to third flux'
                
        mwrfits,values,'observation.fits',hdr,/create

        ;make model fits file

        widget_control,info.t1,get_value=lparam
        pars = lparam(0)
        fixed = lparam(1)
        min = lparam(2)
        max = lparam(3)
        sigma = lparam(4)
        widget_control,info.t2,get_value=sparam
        sdat = sparam
        widget_control,info.t3,get_value=cparam
        cdat = cparam

        sxaddpar,hdr2,'DATE',systime(),'Date of creation'
        sxaddpar,hdr2,'PHI0',pars.phi0,'Luminosity Function Normalization'
        sxaddpar,hdr2,'PHI0_FIX',fixed.phi0,'Fix Phi0 (Y=1/N=0)'
        sxaddpar,hdr2,'PHI0_MIN',min.phi0,'Minimum Phi0 value'
        sxaddpar,hdr2,'PHI0_MAX',max.phi0,'Maximum Phi0 value'
        sxaddpar,hdr2,'PHI0_DP',sigma.phi0,'Sigma Phi0 value'
        sxaddpar,hdr2,'L0',pars.lo,'Luminosity Function Knee'
        sxaddpar,hdr2,'L0_FIX',fixed.lo,'Fix L0 (Y=1/N=0)'
        sxaddpar,hdr2,'L0_MIN',min.lo,'Minimum L0 value'
        sxaddpar,hdr2,'L0_MAX',max.lo,'Maximum L0 value'
        sxaddpar,hdr2,'L0_DP',sigma.lo,'Sigma L0 value'
        sxaddpar,hdr2,'ALPHA',pars.alpha,'Luminosity Function upper slope'
        sxaddpar,hdr2,'ALPHA_FIX',fixed.alpha,'Fix Alpha (Y=1/N=0)'
        sxaddpar,hdr2,'ALPHA_MIN',min.alpha,'Minimum Alpha value'
        sxaddpar,hdr2,'ALPHA_MAX',max.alpha,'Maximum Alpha value'
        sxaddpar,hdr2,'ALPHA_DP',sigma.alpha,'Sigma Alpha value'
        sxaddpar,hdr2,'BETA',pars.beta,'Luminosity Function lower slope'
        sxaddpar,hdr2,'BETA_FIX',fixed.beta,'Fix Beta (Y=1/N=0)'
        sxaddpar,hdr2,'BETA_MIN',min.beta,'Minimum Beta value'
        sxaddpar,hdr2,'BETA_MAX',max.beta,'Maximum Beta value'
        sxaddpar,hdr2,'BETA_DP',sigma.beta,'Sigma Beta value'
        sxaddpar,hdr2,'P',pars.p,'Luminosity Function PHI evolution term'
        sxaddpar,hdr2,'P_FIX',fixed.p,'Fix P (Y=1/N=0)'
        sxaddpar,hdr2,'P_MIN',min.p,'Minimum P value'
        sxaddpar,hdr2,'P_MAX',max.p,'Maximum P value'
        sxaddpar,hdr2,'P_DP',sigma.p,'Sigma P value'
        sxaddpar,hdr2,'Q',pars.q,'Luminosity Function L evolution term'
        sxaddpar,hdr2,'Q_FIX',fixed.q,'Fix Q (Y=1/N=0)'
        sxaddpar,hdr2,'Q_MIN',min.q,'Minimum Q value'
        sxaddpar,hdr2,'Q_MAX',max.q,'Maximum Q value'
        sxaddpar,hdr2,'Q_DP',sigma.q,'Sigma Q value'
        sxaddpar,hdr2,'CEXP',cdat.a0,'Intrinsic luminosity evolution term'
        sxaddpar,hdr2,'CEXP_FIX',cdat.fixed,'Fix CEXP (Y=1/N=0)'
        sxaddpar,hdr2,'CEXP_MIN',cdat.amin,'Minimum CEXP value'
        sxaddpar,hdr2,'CEXP_MAX',cdat.amax,'Maximum CEXP value'
        sxaddpar,hdr2,'CEXP_DP',cdat.sigma,'Sigma CEXP value'
        sxaddpar,hdr2,'RUNS',sdat.runs,'Number of Runs'
        sxaddpar,hdr2,'ZMIN',sdat.zmin,'Minimum Redshift Value'
        sxaddpar,hdr2,'ZMAX',sdat.zmax,'Maximum Redshift Value'
        sxaddpar,hdr2,'DZ',sdat.dz,'Redshit Bin Width'
        sxaddpar,hdr2,'AREA',sdat.area,'Observed Solid Angle'

        templates = [0]
        mwrfits,templates,'model.fits',hdr2,/create

;===============================================================
;Run the actual simulation
;---------------------------------------------------------------
        args = 'observation.fits model.fits sf_templates.fits'
        spawn,'fitter '+args

        read_output
        graphs
     end
     'replot': begin
        read_output
        graphs
     end
     'quit': begin
        widget_control,ev.top,/destroy
     end
     'ot': widget_control,info.ot,get_value=bands
     't1': begin
        widget_control,info.t1,get_value=ldata
        
        i = ev.y
        j = ev.x        

        if (i eq 1) then begin
           if(ldata[i].(j) lt 1) then begin
              ldata[i].(j)=0
              widget_control,ev.id,set_value=ldata
           endif else if(ldata[i].(j) gt 1) then begin
              ldata[i].(j)=1
              widget_control,ev.id,set_value=ldata
           endif
        endif

     end
     't2': widget_control,info.t2,get_value=sdat
     't3': widget_control,info.t3,get_value=cdat
     'info': fit_info
     ELSE:
  ENDCASE

END

pro fit_info
  
  main = widget_base(title="Information about Fitting Methodology",/column)
  t1 = widget_text(main,/wrap,value=['Luminosity Function Fitter Version 1.0' $
  ,'Noah Kurinsky and Anna Sajina' $
  ,'Department of Physics and Astronomy, Tufts University, Medford MA' $
  ,'Last Updated 5/1/13' $
  ,'' $
  ,'Notes on using this interface:' $
  ,'- To allow a parameter to be fitted, set the value of the "fixed" row for that parameter to "0"' $
  ,'- The ranges for the luminosity function parameters are only important for non-fixed parameters' $
  ,'- The "sigma" row refers to the proposed distribution width (stdev of gaussian monte carlo steps'],xsize=80,ysize=20)
  bthold = widget_base(main,/row,/align_center)
  bt = widget_button(bthold,uvalue='close',value='Close',xsize=50,ysize=25)

  widget_control,main,/realize
  xmanager,'fit_info',main,/no_block

end

pro fit_info_event,ev

  widget_control,ev.id,get_uvalue=uvalue

  case uvalue of
     'close': begin
        widget_control,ev.top,/destroy
     end
  endcase

end

pro read_output
  
  !p.thick=5
  !x.thick=5
  !y.thick=5
  !p.charthick=5
  !p.charsize=1.5

  comp = mrdfits('output.fits',0,head,/silent)
  model = mrdfits('output.fits',1,/silent)
  obs = mrdfits('output.fits',2,/silent)

  xysize = fxpar(head,'DIM')
  hist_min = fxpar(head,'H_MIN')
  hist_max = fxpar(head,'H_MAX')
  binsize = fxpar(head,'BINSIZE')
  
  alpha = fxpar(head,'ALPHA')
  beta = fxpar(head,'BETA')
  phi0 = fxpar(head,'PHI0')
  L0 = fxpar(head,'L0')
  p = fxpar(head,'P')
  q = fxpar(head,'Q')

  a = findgen(xysize+1)
  a *= binsize
  a += hist_min

  set_plot,'ps'
  device,filename='model_color_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  
  plot,[hist_min,hist_max],[hist_min,hist_max],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}')
  hb = 0.5

  color = 260*model/(max(model)*1.2)+30
  color[where(model eq 0)] = 0

  loadct,39,/silent

  for i=0,xysize-1 do begin
     for j=0,xysize-1 do begin
        
        xfill = [a[i],a[i],a[i+1],a[i+1]]
        yfill = [a[j],a[j+1],a[j+1],a[j]]

        if(model(i,j) gt 0) then begin
           polyfill,xfill,yfill,color=color(i,j)
        endif

        if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
     endfor
     oplot,[a[i+1],a[i+1]],[hist_min,hist_max],linestyle=1
  endfor

  device,/close

  device,filename='obs_color_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  
  plot,[hist_min,hist_max],[hist_min,hist_max],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}')
  hb = 0.5

  color = 260*obs/(max(obs)*1.2)+30
  color[where(obs eq 0)] = 0

  loadct,39,/silent

  for i=0,xysize-1 do begin
     for j=0,xysize-1 do begin
        
        xfill = [a[i],a[i],a[i+1],a[i+1]]
        yfill = [a[j],a[j+1],a[j+1],a[j]]

        if(obs(i,j) gt 0) then begin
           polyfill,xfill,yfill,color=color(i,j)
        endif

        if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
     endfor
     oplot,[a[i+1],a[i+1]],[hist_min,hist_max],linestyle=1
  endfor

  device,/close

  device,filename='comp_color_hist.eps',xsize=10,ysize=8,/inches,/encapsulated,/times,set_font='Times-Roman',/color
  
  pos = where(comp gt 0)
  neg = where(comp lt 0)
  resmod = comp
  resobs = abs(comp)
  if(neg[0] gt -1) then resmod[neg] = 0
  if(pos[0] gt -1) then resobs[pos] = 0

  plot,[hist_min,hist_max],[hist_min,hist_max],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}')
  hb = 0.5

  for resi=0,1 do begin
     
     case resi of
        0:begin
           wres = resmod
           loadct,3,/silent
        end
        1:begin
           wres = resobs
           loadct,8,/silent
        end
     endcase
     
     color = 249-200*wres/(max(wres))

     for i=0,xysize-1 do begin
        for j=0,xysize-1 do begin
           
           xfill = [a[i],a[i],a[i+1],a[i+1]]
           yfill = [a[j],a[j+1],a[j+1],a[j]]
           
           if(wres(i,j) gt 0) then begin
              polyfill,xfill,yfill,color=color(i,j)
           endif
           
           if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
        endfor
        oplot,[a[i+1],a[i+1]],[hist_min,hist_max],linestyle=1
     endfor
  endfor

  device,/close

  omega_m = 0.28
  omega_l = 0.72
  
  device,filename='sim_lumfunct.eps',xsize=12,ysize=9,/inches,/times,/color,/encapsulated

  plot,[8,13],[1e-10,1e0],/ylog,/nodata,xtitle=textoidl('log_{10}(L_{fir}) [L_{sun}]'),ytitle=textoidl('log_{10}(N(L_{fir})/\Omega dV_c)'),ystyle=1
  lums = indgen(21)/4.0+8.0
  print,lums

  dz = 0.25

  loadct,39,/silent
  for z=dz,5.0,dz do begin
     if z le 2.0 then begin
        t1 = (10^phi0)*((1+z)^p)
        t2 = (10^L0)*((1+z)^q)
     endif else begin
        t1 = (10^phi0)*((3.0)^p)
        t2 = (10^L0)*((3.0)^q)
     endelse

     r = 10^lums/t2
     nsrcs=t1/(r^alpha+r^beta)

     ;ez = sqrt(omega_m*(1+z)^3+omega_l)
     ;dvdz = 4062*lumdist(z)^2*(!pi/180)^2/((1+z)^2*ez)
     ;vol = dvdz*dz
     ;nsrcs*=vol

     ;print,vol
     oplot,lums,nsrcs,color=z*40+20
  endfor

  for i=0.0,5.0,0.1 do begin
     ind = (9+2*i/5)
     oplot,[ind,ind+0.02,ind,ind+0.02,ind,ind+0.02],[10^(-8),10^(-8),10^(-7.9),10^(-7.9),10^(-7.8),10^(-7.8)],psym=2,color=i*40+20,symsize=0.75
  endfor

  loadct,0,/silent
  oplot,[8,13],[1,1],linestyle=1

  xyouts,9.0,10^(-7.2),textoidl('Redshift')
  xyouts,9.0,10^(-7.6),textoidl('0')
  xyouts,10.9,10^(-7.6),textoidl('5')

  device,/close

end

pro graphs

  COMMON simulation_com,info,ldata,ldat0,sdat,cdat,bands

  !p.thick=0
  !x.thick=0
  !y.thick=0
  !p.charthick=0
  !p.charsize=0

  size_screen=get_screen_size()
  size_screen_alt = size_screen*0.85
  size_screen = size_screen*0.8
  gmain = widget_base(title='Simulation Output',/column,xsize=size_screen[0],ysize=size_screen_alt[1])
  r1 = widget_base(gmain,/row) ;,xsize=size_screen[0],ysize=size_screen[1])
  r2 = widget_base(gmain,/row) ;,xsize=size_screen[0],ysize=size_screen[1])
  r2b = widget_base(gmain,/row) ;,xsize=size_screen[0],ysize=size_screen[1])
  r3 = widget_base(gmain,/row) ;,xsize=size_screen[0],ysize=size_screen[1])

  widget_control,gmain,set_uvalue=mnum
  xdim = fix(size_screen[0]/3.0)
  ydim = fix(size_screen[1]/3.0)

  loadct,0,/silent

  lumfunct = widget_draw(r1,xsize=xdim,ysize=ydim)
  redshift = widget_draw(r1,xsize=xdim,ysize=ydim)
  models = widget_draw(r1,xsize=xdim,ysize=ydim)
  dcount1 = widget_draw(r2,xsize=xdim,ysize=ydim)
  dcount2 = widget_draw(r2,xsize=xdim,ysize=ydim)
  dcount3 = widget_draw(r2,xsize=xdim,ysize=ydim)
  sim_colors = widget_draw(r2b,xsize=xdim,ysize=ydim)
  obs_colors = widget_draw(r2b,xsize=xdim,ysize=ydim)
  comp_colors = widget_draw(r2b,xsize=xdim,ysize=ydim)
  
  refresh = widget_button(r3,uvalue='refresh',value='Refresh')
  close = widget_button(r3,uvalue='close',value='Close')
  quit = widget_button(r3,uvalue='quit',value='Quit')
  
  widget_control,gmain,/realize
  xmanager,'graphs',gmain,/no_block
  
  set_plot,'x'
  loadct,0,/silent

  res = mrdfits('output.fits',0,head,/silent)

  alpha = fxpar(head,'ALPHA')
  beta = fxpar(head,'BETA')
  phi0 = fxpar(head,'PHI0')
  L0 = fxpar(head,'L0')
  p = fxpar(head,'P')
  q = fxpar(head,'Q')

  omega_m = 0.28
  omega_l = 0.72
  
  dists = mrdfits('output.fits',3,head,/silent)

  gpts = where(dists.f3 gt 0)
  f1 = dists[gpts].f1
  f2 = dists[gpts].f2
  f3 = dists[gpts].f3
  z = dists[gpts].z
  m = dists[gpts].m
  lum = dists[gpts].lum

  pnum_out = fxpar(head,'tfields')

  ;; for i=6,pnum_out-1 do begin
  ;;    if (i eq 6) then begin
  ;;       hists = histogram(dists[gpts].(i),binsize=0.1,locations=xh,min=0,max=10)
  ;;       hmax = max(hists)
  ;;       xmax = max(xh)
  ;;    endif else begin
  ;;       h = histogram(dists[gpts].(i),binsize=0.1,locations=xh,min=0)
  ;;       if(max(xh) gt xmax) then begin
  ;;          xmax = max(xh)
  ;;          for j=6,i-1 do begin
  ;;             if(i eq 6) then begin
  ;;                hists = histogram(dists[gpts].(i),binsize=0.1,locations=xh,min=0,max=xmax)
  ;;             endif else begin
  ;;                htemp = histogram(dists[gpts].(i),binsize=0.1,locations=xh,min=0,max=xmax)
  ;;                hists = [hists,htemp]
  ;;             endelse
  ;;          endfor
  ;;       endif else begin
  ;;          h = histogram(dists[gpts].(i),binsize=0.1,locations=xh,min=0,max=xmax)
  ;;       endelse
  ;;       temp = max(h)
  ;;       hmax = (temp gt hmax) ? temp : hmax
  ;;       hists = [[hists],[h]]
  ;;    endelse
  ;; endfor     

  widget_control,lumfunct,get_value=index
  wset,index
  device,decomposed=0
  ;h = histogram(lum,nbins=50,locations=xh)
  ;plot,xh,h,psym=10,xstyle=1,yrange=[0.1,1000],ystyle=1,xtitle='Log(Luminosity (W/Hz))',ytitle='dN/(dL/dHz)',/ylog,title='Luminosity Function'

  plot,[8,13],[1e-10,1e0],/ylog,/nodata,xtitle=textoidl('log_{10}(L_{fir}) [L_{sun}]'),ytitle=textoidl('log_{10}(N(L_{fir})/\Omega dV_C)'),ystyle=1
  lums = indgen(21)/4.0+8.0
  print,lums

  dz = 0.5

  loadct,39,/silent
  for zi=dz,5.0,dz do begin
     if zi le 2.0 then begin
        t1 = (10^phi0)*((1+zi)^p)
        t2 = (10^L0)*((1+zi)^q)
     endif else begin
        t1 = (10^phi0)*((3.0)^p)
        t2 = (10^L0)*((3.0)^q)
     endelse

     r = 10^lums/t2
     nsrcs=t1/(r^alpha+r^beta)

     ;ez = sqrt(omega_m*(1+zi)^3+omega_l)
     ;dvdz = 4062*lumdist(zi)^2*(!pi/180)^2/((1+zi)^2*ez)
     ;vol = dvdz*dz
     ;nsrcs*=vol

     ;print,vol
     oplot,lums,nsrcs,color=zi*40+20
  endfor

  for i=0.0,5.0,0.1 do begin
     ind = (9+2*i/5)
     oplot,[ind,ind+0.02,ind,ind+0.02,ind,ind+0.02],[10^(-8),10^(-8),10^(-7.9),10^(-7.9),10^(-7.8),10^(-7.8)],psym=2,color=i*40+20,symsize=0.75
  endfor
  
  loadct,0,/silent
  oplot,[8,13],[1,1],linestyle=1

  xyouts,9.0,10^(-7.2),'Redshift'
  xyouts,9.0,10^(-7.6),'0'
  xyouts,10.9,10^(-7.6),'5'

  widget_control,redshift,get_value=index
  wset,index
  h = histogram(z,binsize=0.1,locations=xh,min=0.2,max=5.0)
  plot,xh,h,psym=10,xrange=[0,5],xstyle=1,xtitle='z',ytitle='dN/dz',title='Redshift Distribution'

  widget_control,models,get_value=index
  wset,index

  h = histogram(alog10(f1/1e3),nbins=50,locations=xh,min=-3,max=1)
  pts = where (h le 0)
  h[pts] = 0.01
  plot,xh,h,psym=10,/ylog,xrange=[-2.2,-0.5],yrange=[1e-1,1e3],xstyle=1,ystyle=1

  ;; help,hists
  ;; for i=6,pnum_out-1 do begin
  ;;    if(i eq 6) then begin
  ;;       plot,xh,hists(*,i-6),psym=10,xstyle=1,xtitle='m',ytitle='dN/dm',xrange=[0,xmax],title='Model Distribution',yrange=[0,hmax*1.2]
  ;;    endif else begin
  ;;       oplot,xh,hists(*,i-6),linestyle=(i-6),psym=10
  ;;    endelse
  ;; endfor

  widget_control,dcount1,get_value=index
  wset,index
  
  ;Herschel ATLAS counts at 250,350 and 500 (Clements et al. 2010)
  if( file_test('counts_clements10.dat')) then begin
     readcol,'counts_clements10.dat',skipline=2,numline=16,flux,nbin,corr,int_counts,int_err,diff_counts,diff_err,/silent
     flux=flux/1.d3
     plot,flux,diff_counts,psym=1,symsize=2,xtitle=TeXtoIDL('F_{250}[Jy]'),ytitle=TeXtoIDL('(dN/dS)S^{2.5} [gal ster^{-1} J^{1.5}]'),/xlog,/ylog,title='Band 1 Counts',yrange=[5.d2,1.d5],ystyle=1
     oploterr,flux,diff_counts,diff_err
  endif else begin
     print,'Error: File "counts_clements10.dat" not found'
     plot,[100,500],[5.d2,1.d5],/nodata,psym=1,symsize=2,xtitle=TeXtoIDL('F_{250}[Jy]'),ytitle=TeXtoIDL('(dN/dS)S^{2.5} [gal ster^{-1} J^{1.5}]'),/xlog,/ylog,title='Band 1 Counts',ystyle=1
  endelse

  h = histogram(alog10(f1/1e3),nbins=50,locations=xh,min=-3,max=1)
  ;I think the area covered here is not
  ;quite right needs to be sorted out better, but for now just scale up
  h=h*2.d3
  ;dlogS=dS/S
  df1=(shift(xh,-1)-xh)*10.^(xh)
  dcounts=(h/df1)*f1^(2.5)
  xh=10.^(xh)
  oplot,xh,dcounts,psym=10

  widget_control,dcount2,get_value=index
  wset,index

  h = histogram(alog10(f2/1e3),nbins=50,locations=xh,min=-3,max=0)
  pts = where(h le 0)
  h[pts] = 0.01
  plot,xh,h,psym=10,/ylog,xrange=[-2.2,-0.5],yrange=[1e-1,1e3],ystyle=1,xstyle=1,xtitle='F_{350}[Log(Jy)]',ytitle='dN/dS (Log)',title='Band 2 Counts'

  widget_control,dcount3,get_value=index
  wset,index

  h = histogram(alog10(f3/1e3),nbins=50,locations=xh,min=-3,max=0)
  pts = where(h le 0)
  h[pts] = 0.01
  plot,xh,h,psym=10,/ylog,xrange=[-2.2,-0.5],yrange=[1e-1,1e3],ystyle=1,xstyle=1,xtitle='F_{500}[Log(Jy)]',ytitle='dN/dS (Log)',title='Band 3 Counts'

  comp = mrdfits(info.oname,0,head,/silent)
  model = mrdfits(info.oname,1,/silent)
  obs = mrdfits(info.oname,2,/silent)

  xysize = fxpar(head,'DIM')
  hist_min = fxpar(head,'H_MIN')
  hist_max = fxpar(head,'H_MAX')
  binsize = fxpar(head,'BINSIZE')
  
  a = findgen(xysize+1)
  a *= binsize
  a += hist_min
  
  loadct,39,/silent
  device,decomposed=0
  widget_control,sim_colors,get_value=index
  wset,index
  
  plot,[hist_min,hist_max],[hist_min,hist_max],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}'),title='Model Color Distribution'
  hb = 0.5

  color = 256*model/(max(model)*1.2)+30
  color[where(model eq 0)] = 0

  for i=0,xysize-1 do begin
     for j=0,xysize-1 do begin
        
        xfill = [a[i],a[i],a[i+1],a[i+1]]
        yfill = [a[j],a[j+1],a[j+1],a[j]]

        if(model(i,j) gt 0) then begin
           polyfill,xfill,yfill,color=color(i,j)
        endif

        if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
     endfor
     oplot,[a[i+1],a[i+1]],[hist_min,hist_max],linestyle=1
  endfor

  widget_control,obs_colors,get_value=index
  wset,index
  
  plot,[hist_min,hist_max],[hist_min,hist_max],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}'),title='Observed Color Distribution'
  hb = 0.5

  color = 256*obs/(max(obs)*1.2)+30
  color[where(obs eq 0)] = 0

  for i=0,xysize-1 do begin
     for j=0,xysize-1 do begin
        
        xfill = [a[i],a[i],a[i+1],a[i+1]]
        yfill = [a[j],a[j+1],a[j+1],a[j]]

        if(obs(i,j) gt 0) then begin
           polyfill,xfill,yfill,color=color(i,j)
        endif

        if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
     endfor
     oplot,[a[i+1],a[i+1]],[hist_min,hist_max],linestyle=1
  endfor

  widget_control,comp_colors,get_value=index
  wset,index
    
  pos = where(comp gt 0)
  neg = where(comp lt 0)
  resmod = comp
  resobs = abs(comp)
  resmod[neg] = 0
  resobs[pos] = 0

  plot,[hist_min,hist_max],[hist_min,hist_max],/nodata,xstyle=1,ystyle=1,xminor=1,yminor=1,xtitle=textoidl('\alpha_{250}^{500}'),ytitle=textoidl('\alpha_{350}^{500}'),title='Color Distribution Comparison'
  hb = 0.5

  for resi=0,1 do begin
     
     case resi of
        0:begin
           wres = resmod
           loadct,3,/silent
        end
        1:begin
           wres = resobs
           loadct,8,/silent
        end
     endcase
     
     color = 249-200*wres/(max(wres))

     for i=0,xysize-1 do begin
        for j=0,xysize-1 do begin
           
           xfill = [a[i],a[i],a[i+1],a[i+1]]
           yfill = [a[j],a[j+1],a[j+1],a[j]]
           
           if(wres(i,j) gt 0) then begin
              polyfill,xfill,yfill,color=color(i,j)
           endif
           
           if(i eq xysize-1) then oplot,[hist_min,hist_max],[a[j+1],a[j+1]],linestyle=1        
        endfor
        oplot,[a[i+1],a[i+1]],[hist_min,hist_max],linestyle=1
     endfor
  endfor

end

pro graphs_event,ev
  COMMON simulation_com,info,ldata,ldat0,sdat,cdat,bands

  widget_control,ev.id,get_uvalue=uvalue

  case uvalue of
     'refresh': begin
        widget_control,ev.top,get_uvalue=mnum
        widget_control,ev.top,/destroy
        graphs
     end
     'close': widget_control,ev.top,/destroy
     'quit': begin
        widget_control,ev.top,/destroy
        widget_control,info.base,/destroy
     end
  endcase

end
