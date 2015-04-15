 pro make_template_file

data=LoadData(2)
red=GetColor('Red',5)
blue=GetColor('Blue',6)
;==========================================================================
;start with the default Rieke et al 2009 starforming galaxy IR SED templates
;---------------------------------------------------------------------------
;the lambda units are microns
;the template units are W/Hz (i.e. L_nu) 
;all templates are given in the restframe

rieke_templ=mrdfits('sf_templates.fits')
old_header=headfits('sf_templates.fits')
lambda=rieke_templ[*,0] ;in microns
tmp=size(rieke_templ)
nlam=tmp[1]
ntmpl=tmp[2]

lums=[9.75,10.00,10.25,10.50,10.75,11.00,11.25,11.50,11.75,12.00,12.25,12.50,12.75,13.00]
nlum=ntmpl-1
print,'Number of templates Rieke',nlum
;============================================================================
;check templates
;----------------------------------------------------------------------------

templ1=rieke_templ[*,1]

plot,lambda,templ1,/xlog,/ylog,xtitle='lambda [microns]',ytitle='L_nu [W/Hz]',xrange=[0.1,5.d5],yrange=[1d20,1d28],ystyle=1

oplot,lambda,templ1,color=red

templ10=rieke_templ[*,10]
oplot,lambda,templ10,color=blue

;need to stick a stellar component to these SF templates
readcol,'/Users/annie/code/idl/auxdata/swire/M82_template_norm.sed',lam_m82,flam_m82,format='f,f'

lam_m82=lam_m82/1.d4 ;convert to microns from Angstroms
nonzero=where(templ1 gt 0)
diff=abs(lam_m82-lambda[nonzero[0]])
stitch=where(diff eq min(diff))
print,'lambda stitch at: ',lambda[nonzero[0]],lam_m82[stitch]

fnu_m82=(flam_m82*lam_m82^2.)
;print,templ1[nonzero[0]],fnu_m82[stitch]
scale=templ1[nonzero[0]]/fnu_m82[stitch]
fnu_m82_scaled=scale#fnu_m82
oplot,lam_m82[0:stitch-1],fnu_m82_scaled[0:stitch-1],color=red

nonzero=where(templ10 gt 0)
diff=abs(lam_m82-lambda[nonzero[0]])
stitch=where(diff eq min(diff))
scale=templ10[nonzero[0]]/fnu_m82[stitch]
fnu_m82_scaled=scale#fnu_m82
oplot,lam_m82[0:stitch-1],fnu_m82_scaled[0:stitch-1],color=blue

lambda_nonzero=lambda[nonzero[0]:n_elements(lambda)-1]
lambda_new=[lam_m82[0:stitch-1],lambda_nonzero]
nlam=n_elements(lambda_new)
dlambda=lambda_new-shift(lambda_new,1) ;in microns
dnu=(dlambda*(3.d8/lambda_new^2.))*1.d6 ;to get it to Hz

rieke_templ_new=fltarr(n_elements(lambda_new),ntmpl)

rieke_templ_new[*,0]=lambda_new 
for itempl=1,ntmpl-1 do begin
   tmp_templ=rieke_templ[*,itempl]
   nonzero=where(tmp_templ gt 0)
   diff=abs(lam_m82-lambda[nonzero[0]])
   stitch=where(diff eq min(diff))
   lambda_nonzero=lambda[nonzero[0]:n_elements(lambda)-1]
   lambda_tmp=[lam_m82[0:stitch-1],lambda_nonzero]
   scale=tmp_templ[nonzero[0]]/fnu_m82[stitch]
   fnu_m82_scaled=scale#fnu_m82
   tmp_templ_nonzero=tmp_templ[nonzero[0]:n_elements(tmp_templ)-1]
   tmp_new=[fnu_m82_scaled[0:stitch-1],tmp_templ_nonzero]
   rieke_templ_new[*,itempl]=interpol(tmp_new,lambda_tmp,lambda_new)
endfor

;templ1_new=rieke_templ_new[*,10]
;oplot,lambda_new,rieke_templ_new[*,1],color=red
;oplot,lambda_new,rieke_templ_new[*,10],color=blue

;Rieke et al. 2009, use L_TIR which is 5-1000um
gpts=where((lambda_new ge 5.0) and (lambda_new lt 1000.0))
;lum1=alog10(total(templ1_new[gpts]*dnu[gpts])/3.846d26)
;print,lum1 ;here I am getting luminosities that are close to but not exactly as given for the templates (vary by +/-0.13), not sure why the difference really

;============================================================================
;read-in new templates and cast in appropriate format
;----------------------------------------------------------------------------

useset=1 ;Sajina et al. 2012
useset=2 ;Kirkpatrick et al. 2015

if(useset eq 1) then begin
   y_z1_comp=fltarr(nlam,nlum)
   y_z2_comp=fltarr(nlam,nlum)
   y_z1_lowtau=fltarr(nlam,nlum)
   y_z2_lowtau=fltarr(nlam,nlum)
   y_z1_hightau=fltarr(nlam,nlum)
   y_z2_hightau=fltarr(nlam,nlum)

;agn_lowtau (Sajina et al. 2012)
   readcol,'IRSsupersample_comp_template.dat',lam,y1,y1_low,y1_high,y2,y2_low,y2_high,skipline=8
;plot,lam,y1,/xlog,/ylog,yrange=[0.01,10],ystyle=1,xrange=[0.5,1000],xstyle=1
;oplot,lam*(1+4),y1,linestyle=2,color=red
;just testing, overplot the wise filters
;filter_data=mrdfits('/Users/annie/students/noah_kurinsky/SurveySim/trunk/model/model.fits',1)
;help,filter_data,/struct
;oplot,filter_data.lambda1/1.d4,filter_data.transmission1,linestyle=1
;print,max(filter_data.lambda1),max(filter_data.transmission1)
;oplot,filter_data.lambda2/1.d4,filter_data.transmission2,linestyle=1
;oplot,filter_data.lambda3/1.d4,filter_data.transmission3,linestyle=1

   nu_my=(3.d8/lam)*1.d6        ;in Hz
   y_tmp=interpol(y1/nu_my,lam,lambda_new)
   tmp_lum=alog10(total(y_tmp[gpts]*dnu[gpts])/3.846d26)
   for i=0,nlum-1 do begin
      ltir=lums[i]
      scale=10.^(ltir-tmp_lum)
      y_z1_comp[*,i]=y_tmp*scale
   endfor
   y_tmp=interpol(y2/nu_my,lam,lambda_new)
   tmp_lum=alog10(total(y_tmp[gpts]*dnu[gpts])/3.846d26)
   for i=0,nlum-1 do begin
      ltir=lums[i]
      scale=10.^(ltir-tmp_lum)
      y_z2_comp[*,i]=y_tmp*scale
   endfor

   readcol,'IRSsupersample_agn_template_lowtau.dat',lam,y1,y1_low,y1_high,y2,y2_low,y2_high,skipline=8
   nu_my=(3.d8/lam)*1.d6        ;in Hz
   y_tmp=interpol(y1/nu_my,lam,lambda_new)
   tmp_lum=alog10(total(y_tmp[gpts]*dnu[gpts])/3.846d26)
   for i=0,nlum-1 do begin
      ltir=lums[i]
      scale=10.^(ltir-tmp_lum)
      y_z1_lowtau[*,i]=y_tmp*scale
   endfor
   y_tmp=interpol(y2/nu_my,lam,lambda_new)
   tmp_lum=alog10(total(y_tmp[gpts]*dnu[gpts])/3.846d26)
   for i=0,nlum-1 do begin
      ltir=lums[i]
      scale=10.^(ltir-tmp_lum)
      y_z2_lowtau[*,i]=y_tmp*scale
   endfor

;agn_hightau (Sajina et al. 2012)
   readcol,'IRSsupersample_agn_template_hightau.dat',lam,y1,y1_low,y1_high,y2,y2_low,y2_high,skipline=8
   nu_my=(3.d8/lam)*1.d6        ;in Hz
   y_tmp=interpol(y1/nu_my,lam,lambda_new)
   tmp_lum=alog10(total(y_tmp[gpts]*dnu[gpts])/3.846d26)
   for i=0,nlum-1 do begin
      ltir=lums[i]
      scale=10.^(ltir-tmp_lum)
      y_z1_hightau[*,i]=y_tmp*scale
   endfor
   y_tmp=interpol(y2/nu_my,lam,lambda_new)
   tmp_lum=alog10(total(y_tmp[gpts]*dnu[gpts])/3.846d26)
   for i=0,nlum-1 do begin
      ltir=lums[i]
      scale=10.^(ltir-tmp_lum)
      y_z2_hightau[*,i]=y_tmp*scale
   endfor

   nz=3 ;z=0 (the Rieke templates),z=1 (or rather some appropriate bin), and z=2 (again some bin)
   nsed=4 ;assume at least 4 choices of SED per L-z bin (potentially 2 SF (or 1 SF, 1 composite) and 2 types of AGN (more or less obscured)
   templates=fltarr(nlam,ntmpl,nz,nsed)
   for i=0,ntmpl-1 do begin
      for j=0,nz-1 do begin
                                ;for now, starforming templates constant with redshift
      templates[*,i,j,0]=rieke_templ_new[*,i] ;note this includes the lambda array
      if(i gt 0) then begin
         if(j eq 1) then begin
            templates[*,i,j,1]=y_z1_comp[*,i-1]
            templates[*,i,j,2]=y_z1_lowtau[*,i-1]
            templates[*,i,j,3]=y_z1_hightau[*,i-1]
         endif
         if(j eq 2) then begin
            templates[*,i,j,1]=y_z2_comp[*,i-1]
            templates[*,i,j,2]=y_z2_lowtau[*,i-1]
            templates[*,i,j,3]=y_z2_hightau[*,i-1]
         endif
      endif
   endfor
   endfor

;oplot,lambda_new,templates[*,1,1,1]
;modify the header
   my_header=old_header
;note these are placeholders, need to be sorted out exactly -- plus
;                                                              these
;are ranges rather than fixed points
   SXADDPAR, my_header, 'Z0', 0,'lowest redshift template'
   SXADDPAR, my_header, 'Z1', 1,'middle redshift template'
   SXADDPAR, my_header, 'Z2', 2,'highest redshift template'

   writefits,'Sajinaetal2012_templates.fits',templates,my_header

endif

if(useset eq 2) then begin
   y_z1_comp=fltarr(nlam,nlum)
   y_z2_comp=fltarr(nlam,nlum)
   y_z1_lowtau=fltarr(nlam,nlum)
   y_z2_lowtau=fltarr(nlam,nlum)
   y_z1_hightau=fltarr(nlam,nlum)
   y_z2_hightau=fltarr(nlam,nlum)




   nz=3 ;z=0 (the Rieke templates),z=1 (or rather some appropriate bin), and z=2 (again some bin)
   nsed=4 ;assume at least 4 choices of SED per L-z bin (potentially 2 SF (or 1 SF, 1 composite) and 2 types of AGN (more or less obscured)
   templates=fltarr(nlam,ntmpl,nz,nsed)
   for i=0,ntmpl-1 do begin
      for j=0,nz-1 do begin
                                ;for now, starforming templates constant with redshift
      templates[*,i,j,0]=rieke_templ_new[*,i] ;note this includes the lambda array
      if(i gt 0) then begin
         if(j eq 1) then begin
            templates[*,i,j,1]=y_z1_comp[*,i-1]
            templates[*,i,j,2]=y_z1_lowtau[*,i-1]
            templates[*,i,j,3]=y_z1_hightau[*,i-1]
         endif
         if(j eq 2) then begin
            templates[*,i,j,1]=y_z2_comp[*,i-1]
            templates[*,i,j,2]=y_z2_lowtau[*,i-1]
            templates[*,i,j,3]=y_z2_hightau[*,i-1]
         endif
      endif
   endfor
   endfor

;oplot,lambda_new,templates[*,1,1,1]
;modify the header
   my_header=old_header
;note these are placeholders, need to be sorted out exactly -- plus
;                                                              these
;are ranges rather than fixed points
   SXADDPAR, my_header, 'Z0', 0,'lowest redshift template'
   SXADDPAR, my_header, 'Z1', 1,'middle redshift template'
   SXADDPAR, my_header, 'Z2', 2,'highest redshift template'

   writefits,'Kirkpatricketal2015_templates.fits',templates,my_header
endif

skiptoend:

end
