pro make_template_file

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
;============================================================================
;check templates
;----------------------------------------------------------------------------

templ1=rieke_templ[*,10]

plot,lambda,templ1,/xlog,/ylog,xtitle='lambda [microns]',ytitle='L_nu [W/Hz]'

dlambda=lambda-shift(lambda,1) ;microns
dnu=(dlambda*(3.d8/lambda^2.))*1.d6 ;to get it to Hz

;Rieke et al. 2009, use L_TIR which is 5-1000um
gpts=where((lambda ge 5.0) and (lambda lt 1000.0))
lum1=alog10(total(templ1[gpts]*dnu[gpts])/3.846d26)
print,lum1 ;here I am getting luminosities that are close to but not exactly as given for the templates (vary by +/-0.13), not sure why the difference really

;============================================================================
;read-in new templates and cast in appropriate format
;----------------------------------------------------------------------------
y_z1_comp=fltarr(nlam,nlum)
y_z2_comp=fltarr(nlam,nlum)
y_z1_lowtau=fltarr(nlam,nlum)
y_z2_lowtau=fltarr(nlam,nlum)
y_z1_hightau=fltarr(nlam,nlum)
y_z2_hightau=fltarr(nlam,nlum)

;agn_lowtau (Sajina et al. 2012)
readcol,'IRSsupersample_comp_template.dat',lam,y1,y1_low,y1_high,y2,y2_low,y2_high,skipline=8
nu_my=(3.d8/lam)*1.d6 ;in Hz
y_tmp=interpol(y1/nu_my,lam,lambda)
tmp_lum=alog10(total(y_tmp[gpts]*dnu[gpts])/3.846d26)
for i=0,nlum-1 do begin
   ltir=lums[i]
   scale=10.^(ltir-tmp_lum)
   y_z1_comp[*,i]=y_tmp*scale
endfor
y_tmp=interpol(y2/nu_my,lam,lambda)
tmp_lum=alog10(total(y_tmp[gpts]*dnu[gpts])/3.846d26)
for i=0,nlum-1 do begin
   ltir=lums[i]
   scale=10.^(ltir-tmp_lum)
   y_z2_comp[*,i]=y_tmp*scale
endfor

readcol,'IRSsupersample_agn_template_lowtau.dat',lam,y1,y1_low,y1_high,y2,y2_low,y2_high,skipline=8
nu_my=(3.d8/lam)*1.d6 ;in Hz
y_tmp=interpol(y1/nu_my,lam,lambda)
tmp_lum=alog10(total(y_tmp[gpts]*dnu[gpts])/3.846d26)
for i=0,nlum-1 do begin
   ltir=lums[i]
   scale=10.^(ltir-tmp_lum)
   y_z1_lowtau[*,i]=y_tmp*scale
endfor
y_tmp=interpol(y2/nu_my,lam,lambda)
tmp_lum=alog10(total(y_tmp[gpts]*dnu[gpts])/3.846d26)
for i=0,nlum-1 do begin
   ltir=lums[i]
   scale=10.^(ltir-tmp_lum)
   y_z2_lowtau[*,i]=y_tmp*scale
endfor

;agn_hightau (Sajina et al. 2012)
readcol,'IRSsupersample_agn_template_hightau.dat',lam,y1,y1_low,y1_high,y2,y2_low,y2_high,skipline=8
nu_my=(3.d8/lam)*1.d6 ;in Hz
y_tmp=interpol(y1/nu_my,lam,lambda)
tmp_lum=alog10(total(y_tmp[gpts]*dnu[gpts])/3.846d26)
for i=0,nlum-1 do begin
   ltir=lums[i]
   scale=10.^(ltir-tmp_lum)
   y_z1_hightau[*,i]=y_tmp*scale
endfor
y_tmp=interpol(y2/nu_my,lam,lambda)
tmp_lum=alog10(total(y_tmp[gpts]*dnu[gpts])/3.846d26)
for i=0,nlum-1 do begin
   ltir=lums[i]
   scale=10.^(ltir-tmp_lum)
   y_z2_hightau[*,i]=y_tmp*scale
endfor

;============================================================================
;make new templates.fits file
;----------------------------------------------------------------------------

nz=3 ;z=0 (the Rieke templates),z=1 (or rather some appropriate bin), and z=2 (again some bin)
nsed=4 ;assume at least 4 choices of SED per L-z bin (potentially 2 SF (or 1 SF, 1 composite) and 2 types of AGN (more or less obscured)
templates=fltarr(nlam,ntmpl,nz,nsed)
for i=0,ntmpl-1 do begin
   for j=0,nz-1 do begin
      ;for now, starforming templates constant with redshift
      templates[*,i,j,0]=rieke_templ[*,i] ;note this includes the lambda array
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

;modify the header
my_header=old_header
writefits,'my_templates.fits',templates,my_header
skiptoend:

end
