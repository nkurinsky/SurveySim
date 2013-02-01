PRO runmodels_event, ev
COMMON database, info, frame,lum_0,phi_0,dphi_0,zs,lcube,sedpar,lam,ext,powpah,filts,iter

widget_control, /draw_button_events, get_uvalue=type, ev.id
thisEvent = Tag_Names(ev, /Structure_Name)

IF thisEvent EQ "WIDGET_TRACKING" THEN $
  widget_control, ev.id, /INPUT_FOCUS

IF thisEvent NE "WIDGET_TRACKING" THEN BEGIN
    CASE type OF
        "GO": BEGIN
            plot_ccd
        END
        "QUIT": BEGIN
            widget_control, /reset
        END
        ELSE:                   ;print,'unknown event'
    ENDCASE
ENDIF

END


PRO runmodels
COMMON database, info, frame,lum_0,phi_0,dphi_0,zs,lcube,sedpar,lam,ext,powpah,filts,iter

iter=1

lf_local = read_ascii('serjeant_lf.dat')
lum_0 = lf_local.field1[0,*]+0.5
phi_0 = lf_local.field1[3,*]*10.^(lf_local.field1[4,*])
dphi_0 = lf_local.field1[5,*]*10.^(lf_local.field1[6,*])

zs=fltarr(20)
for iz=0,19 do zs[iz]=0.0+iz*0.2
lcube=fltarr(20,n_elements(lum_0),2)

info = {$
       base:0L, $
       base1:0L, $
       base2:0L, $
       base3:0L, $
       base4:0L, $
       base5:0L, $
;       spectrum:0L, $
       text:strarr(3), $
       win_id1:0L, $
       win_id2:0L, $
       win_id3:0L, $
       win_id4:0L, $
       win_id5:0L, $
       win_id6:0L, $
       cursor:[0.,0.], $
       xsize1:128., $
       ysize1:128., $
       xsize2:128., $
       ysize2:128., $
       magnification:2., $
       minscale:-2, $
       maxscale:2,  $ 
       draw1:0L, $
       draw2:0L, $
       draw3:0L, $
       draw4:0L, $
       draw5:0L, $
       draw6:0L, $
       table:0L, $
       ncolors:0L,$
       tmp1:0.,$ ;; temporal variable
       tmp2:0. }

;initialize SED parameters
sedpar = {$
;       z:2.0, $
       agnfrac:0.5, $
       lbol:12.8, $
       mstar:11.2, $
       tauv:0.0, $
       alpha_hot:1.8, $
       alpha:4.0, $
       tau2:15.0, $
       pfrac:0.6, $  
       td: 50.0 $ 
}

restore,'intervectors.save'
lam=lam_new
ext=ext_new ;mw draine&li at this point-may include different things here?
powpah=powpah_new

inps=read_ascii('inputs_hers.dat')
lam_tmp=10.^(inps.field01[0,*])
filt24_mips=inps.field01[14,*]
filt70_mips=inps.field01[19,*]
filt100_pacs=inps.field01[20,*]
filt250_spire=inps.field01[23,*]
filt350_spire=inps.field01[24,*]
filt500_spire=inps.field01[26,*]

filt24=interpol(filt24_mips,lam_tmp,lam)
filt100=interpol(filt100_pacs,lam_tmp,lam)
filt250=interpol(filt250_spire,lam_tmp,lam)
filt350=interpol(filt350_spire,lam_tmp,lam)
filt500=interpol(filt500_spire,lam_tmp,lam)

filts= { $
         r:filtr, $
         i:filti, $
         irac1:filt3, $
         irac2:filt4, $
         irac3:filt5, $
         irac4:filt8, $
         mips1:filt24, $
         mips2:filt70, $
         pacs100:filt100, $
         s170:filt170, $
         spire250:filt250, $
         spire350:filt350, $
         spire500:filt500, $
         s850:filt850 }

LaunchWidgets

lfbox
seds
seltable

END

PRO seltable_event, ev  
COMMON database,info,frame,lum_0,phi_0,dphi_0,zs,lcube,sedpar,lam,ext,powpah,filts,iter

widget_control,info.table,get_value=partable,use_table_select=[0,1,7,1]
;sedpar.z=partable[0]
;sedpar.lbol=partable[1]
;sedpar.mstar=partable[2]
;sedpar.tauv=partable[3]
;sedpar.alpha=partable[4]
;sedpar.tau2=partable[5]
;sedpar.pfrac=partable[6]
;sedpar.td=partable[7]

;sedpar.z=partable[0]
sedpar.agnfrac=partable[0]
sedpar.lbol=partable[1]
sedpar.mstar=partable[2]
;sedpar.tauv=partable[3]
sedpar.alpha_hot=partable[3]
sedpar.alpha=partable[4]
sedpar.tau2=partable[5]
sedpar.pfrac=partable[6]
sedpar.td=partable[7]

iter=iter+1
seds

END

PRO seltable_quit_event, ev  
WIDGET_CONTROL, ev.TOP, /DESTROY  
END
  
PRO seltable
COMMON database,info,frame,lum_0,phi_0,dphi_0,zs,lcube,sedpar,lam,ext,powpah,filts,iter

;d0={param:'z',value:sedpar.z}
;d1={param:'Lir',value:sedpar.lbol}
;d2={param:'Mstar',value:sedpar.mstar}
;d3={param:'tau_*',value:sedpar.tauv}
;d4={param:'alpha',value:sedpar.alpha}
;d5={param:'tau_IR',value:sedpar.tau2}
;d6={param:'pfrac',value:sedpar.pfrac}
;d7={param:'Td',value:sedpar.td}

d0={param:'AGNfrac',value:sedpar.agnfrac}
d1={param:'Lir',value:sedpar.lbol}
d2={param:'Mstar',value:sedpar.mstar}
d3={param:'alpha_hot',value:sedpar.alpha_hot}
;d3={param:'tau_*',value:sedpar.tauv}
d4={param:'alpha',value:sedpar.alpha}
d5={param:'tau_V',value:sedpar.tau2}
d6={param:'pfrac',value:sedpar.pfrac}
d7={param:'Td',value:sedpar.td}

data=[d0,d1,d2,d3,d4,d5,d6,d7]

labels = ['Parameter', 'Value']  
max_strlen = strlen('Orbit Radius (AU)')  
maxwidth = max_strlen * !d.x_ch_size + 6   ; ... + 6 for padding  

base = WIDGET_BASE(/COLUMN)  
info.table = WIDGET_TABLE(base, VALUE=data, /COLUMN_MAJOR, $  
ROW_LABELS=labels, COLUMN_LABELS='', $  
COLUMN_WIDTHS=maxwidth, /RESIZEABLE_COLUMNS)  

b_quit = WIDGET_BUTTON(base, VALUE='Close Table', $ 
EVENT_PRO='seltable_quit_event')    
WIDGET_CONTROL, base,/REALIZE   
col_widths = WIDGET_INFO(info.table, /COLUMN_WIDTHS)  

WIDGET_CONTROL, info.table, COLUMN_WIDTHS=col_widths[0], $  
USE_TABLE_SELECT=[-1,-1,3,3]  

widget_control,info.table,/editable,use_table_select=[0,1,6,1]

XMANAGER, 'seltable', base  


END

PRO LaunchWidgets
COMMON database,info,frame,lum_0,phi_0,dphi_0,zs,lcube,sedpar,lam,ext,powpah,filts,iter

;; Colors
device, pseudo_color=8, decompose=0
loadct,3,/silent
tvlct,red,green,blue,/get
info.ncolors=!d.TABLE_SIZE
red[info.ncolors-1]=0B
green[info.ncolors-1]=255B
blue[info.ncolors-1]=0B
tvlct,red,green,blue

;; Screen size
size_screen = get_screen_size()
size_screen=size_screen*0.8
spectrum_xsize=size_screen(0)*2./3. & spectrum_ysize=size_screen(1)/5.

plotSize = 100.
info.magnification = size_screen[1]/(plotSize+info.ysize1)*0.9

;;; Start the interactive part (widgets)
info.base = widget_base(title='Galaxy evolution model',/row)
info.base1= widget_base(info.base, /column)
info.draw1 = widget_draw(info.base1, xsize=info.magnification*info.xsize1,$
                         ysize=info.magnification*info.ysize1 $
                   ,uvalue="DRAW_WINDOW1",retain=2 $
                   ,/button_events, keyboard_events=1,/tracking_events)
info.draw2 = widget_draw(info.base1, xsize=info.magnification*info.xsize1,$
                         ysize=info.magnification*plotSize $
                         ,uvalue="DRAW_WINDOW2",retain=2 $
                         ,/button_events, keyboard_events=1,/tracking_events)
info.base2= widget_base(info.base, /column)
info.draw3 = widget_draw(info.base2, xsize=info.magnification*info.xsize2,$
                         ysize=info.magnification*info.ysize1 $
                         ,uvalue="DRAW_WINDOW3",retain=2 $
                         ,/button_events, keyboard_events=1,/tracking_events)
info.draw4 = widget_draw(info.base2, xsize=info.magnification*info.xsize2,$ 
                         ysize=info.magnification*plotSize $
                         ,uvalue="DRAW_WINDOW4",retain=2 $
                         ,/button_events, keyboard_events=1,/tracking_events)

info.base3=widget_base(info.base1,/row)
info.base4=widget_base(info.base,/column)
info.draw5 = widget_draw(info.base4, xsize=info.magnification*info.xsize2,$
                         ysize=info.magnification*info.ysize1 $
                         ,uvalue="DRAW_WINDOW5",retain=2 $
                         ,/button_events, keyboard_events=1,/tracking_events)
info.draw6 = widget_draw(info.base4, xsize=info.magnification*info.xsize2,$ 
                         ysize=info.magnification*plotSize $
                         ,uvalue="DRAW_WINDOW6",retain=2 $
                         ,/button_events, keyboard_events=1,/tracking_events)

OK = WIDGET_BUTTON(info.base3,UVALUE='GO',VALUE='GO!',xsize=60,ysize=60)
QUIT = WIDGET_BUTTON(info.base3,UVALUE='QUIT',VALUE='QUIT',xsize=60,ysize=60)

widget_control, /realize, info.base,xoffset=0, yoffset=0
widget_control, info.draw1, get_value=win_id1 & info.win_id1=win_id1
widget_control, info.draw2, get_value=win_id2 & info.win_id2=win_id2
widget_control, info.draw3, get_value=win_id3 & info.win_id3=win_id3
widget_control, info.draw4, get_value=win_id4 & info.win_id4=win_id4
widget_control, info.draw5, get_value=win_id5 & info.win_id5=win_id5
widget_control, info.draw6, get_value=win_id6 & info.win_id6=win_id6

xmanager, 'runmodels', info.base, /NO_BLOCK

END

PRO lfbox
COMMON database,info,frame,lum_0,phi_0,dphi_0,zs,lcube,sedpar,lam,ext,powpah,filts,iter

wset, info.win_id1

data=LoadData(2)
white =GetColor('White',1)
black = GetColor('Black',2)
red = GetColor('Red',3)
blue = GetColor('Blue',4)

lum_1 = [11.0,11.6,12.0,12.6]
phi_1 = [-2.52,-2.70,-3.52,-4.30]
phi_1 = 10.^phi_1

lum_2 = [12.0,12.6,13.00,13.6]
phi_2 = [-3.398,-3.398,-4.00,-5.00]
phi_2 = 10.^phi_2

plot,lum_0,phi_0,/nodata,/ylog,xrange=[10,13.7],yrange=[10.^(-7.5),10.^(-2)],xtitle=TeXtoIDL('log(L_{IR})'),ytitle=TeXtoIDL('\phi [Mpc^{-3} dlogL^{-1}]'),color=white,charsize=1.5,title='INPUT'

;Serjeant et al. data
;oploterr,lum_0,phi_0,dphi_0
;oploterr,[10.2],[10.^(-5)],[2.0*10.^(-6)]
;xyouts,[10.3],[10.^(-5)],TeXtoIDL('z=0 (Serjeant et al.)')

z=0.4
lumz = [9.0, 9.5, 10.0, 10.5]
lumz = lumz+1.1
phiz = [5.8*10.^(-3), 2.93*10.^(-3), 5.44*10.^(-4), 8.37*10.^(-5)]
dphiz = [3.25*10.^(-3), 1.65*10.^(-3), 2.37*10.^(-4), 11.66*10.^(-5)]

;oplot,lumz,phiz,color=red

;z=0.5
lumz = [9.5, 10.0, 10.5, 11.0]
lumz = lumz+1.1
phiz = [2.62*10.^(-3), 1.75*10.^(-3), 3.78*10.^(-4), 3.44*10.^(-5)]
dphiz = [1.43*10.^(-3), 0.93*10.^(-3), 2.18*10.^(-4), 4.1*10.^(-5)]

;oplot,lumz,phiz,color=red

;z=0.7
lumz = [10.0, 10.5, 11.0, 11.5]
lumz = lumz+1.1
phiz = [2.9*10.^(-3), 7.72*10.^(-4), 3.35*10.^(-5), 8.38*10.^(-6)]
dphiz = [1.43*10.^(-3), 0.93*10.^(-4), 2.18*10.^(-5), 4.1*10.^(-6)]

;oplot,lumz,phiz,color=red

;z=0.9
lumz =[10.0, 10.5, 11.0, 11.5]
lumz = lumz+1.2
phiz = [2.27*10.^(-3), 1.4*10.^(-3), 1.06*10.^(-4), 1.07*10.^(-5)]
dphiz = [1.43*10.^(-3), 0.93*10.^(-3), 2.18*10.^(-4), 4.1*10.^(-5)]

;oplot,lumz,phiz,color=red

;z=1.1
lumz = [10.0, 10.5, 11.0, 11.5]
lumz = lumz+1.1
phiz = [2.38*10.^(-3), 1.21*10.^(-3), 2.54*10.^(-4), 2.35*10.^(-5)]
dphiz = [1.43*10.^(-3), 0.93*10.^(-3), 2.18*10.^(-4), 4.1*10.^(-5)]

;oplot,lumz,phiz,color=red

;z=1.0 (Caputi et al. 2007)
lumz=[11.5,11.75,12.0,12.2,12.4]
lphiz=[-2.5,-2.7,-3.0,-3.3,-3.9]
phiz=10.^(lphiz)

;oplot,lumz,phiz,color=red,linestyle=1

;z=2.0 (Caputi et al.)
;lumz = [11.5,12.0,12.0,12.5,12.7]
;lphiz = [-2.9,-3.2,-3.4,-3.8,-4.5]
;the first point is stacked
lumz=[11.5,12.1,12.25,12.4,12.7]
lphiz=[-2.9,-3.2,-3.45,-3.8,-4.4]
phiz=10.^(lphiz)

;oplot,lumz,phiz,color=red,linestyle=2

;AGN from our LF paper
;lum8=[11.4,11.8,12.2,12.6]
;lumz=2.79+0.83*lum8 ;upper limit only!
;lphiz=[-4.33,-5.01,-5.42,-6.99]
;phiz=10.^(lphiz)
;oplot,lumz,phiz,psym=sym(1)

;oplot,[10.2,10.5],[5.5*10.^(-7),5.5*10.^(-7)],color=red,linestyle=1
;xyouts,[10.6],[5.*10.^(-7)],TeXtoIDL('z=1 (Caputi et al.)'),color=red
;oplot,[10.2,10.5],[2.5*10.^(-7),2.5*10.^(-7)],color=red,linestyle=2
;xyouts,[10.6],[2.*10.^(-7)],TeXtoIDL('z=2 (Caputi et al.)'),color=red

oplot,[10.2,10.5],[10.5*10.^(-7),10.5*10.^(-7)],color=white,linestyle=0
xyouts,[10.6],[10.5*10.^(-7)],TeXtoIDL('z=0'),color=white
oplot,[10.2,10.5],[5.5*10.^(-7),5.5*10.^(-7)],color=white,linestyle=1
xyouts,[10.6],[5.*10.^(-7)],TeXtoIDL('z=1'),color=white
oplot,[10.2,10.5],[2.5*10.^(-7),2.5*10.^(-7)],color=white,linestyle=2
xyouts,[10.6],[2.5*10.^(-7)],TeXtoIDL('z=2'),color=white
oplot,[10.2,10.5],[1.0*10.^(-7),1.0*10.^(-7)],color=white,linestyle=3
xyouts,[10.6],[1.0*10.^(-7)],TeXtoIDL('z=3'),color=white


;****************************************************************************

alpha1=0.2
alpha2 = -2.2
lstar = 10.8
pstar= -2.5

phi_sm = pstar-alpha1*lum_0+alpha1*lstar-((10.^lum_0)/10.^lstar)*alog10(2.718)

plregime = where(lum_0 gt 11.0)
phi_sm[plregime] = pstar-0.4-alpha2*lstar+alpha2*lum_0[plregime]

phi_sm=10.^(phi_sm)

for iz=0,19 do begin
z = zs[iz]
alpha1=0.2
alpha2 = -2.2 ;+1.3*alog10(1.+z)
lstar_2 = lstar+3.7*alog10(1.+z)
;if(z le 1.0) then lstar_2 = lstar+5.5*alog10(1.+z)
if(z gt 1.0) then lstar_2=lstar+3.7*alog10(1.+1.0)+0.5*alog10(1.+z)
pstar= -2.5
if((z gt 1.0)) then pstar = -2.5-0.8*alog10(1.+z)
if((z gt 2.0)) then pstar=-2.5-1.0*alog10(1.+z)

good = where(lum_0 ge (lstar_2-0.6))
phi_sm_2 = pstar-alpha1*lum_0+alpha1*lstar_2-((10.^lum_0)/10.^lstar_2)*alog10(2.718)
plregime = where(lum_0 gt (lstar_2+0.3))

phi_sm_2[plregime] = pstar-0.4-alpha2*lstar_2+alpha2*lum_0[plregime]

phi_sm_2=10.^(phi_sm_2)

;oplot,lum_0[good],phi_sm_2[good],color=blue
if(z eq 0.0) then oplot,lum_0[good],phi_sm_2[good],color=white
if(z eq 1.0) then oplot,lum_0[good],phi_sm_2[good],color=white,linestyle=1
if(z eq 2.0) then oplot,lum_0[good],phi_sm_2[good],color=white,linestyle=2
if(z eq 3.0) then oplot,lum_0[good],phi_sm_2[good],color=white,linestyle=3
;if(z eq 4.0) then oplot,lum_0[good],phi_sm_2[good],color=white,linestyle=2

;print,z,ageuniv(z)

;dv=dvol(z)
dv=dvol(z,0.1,1.0) ;this assumes per steradian
dlogl=0.125
for il=0,n_elements(lum_0)-1 do begin
;SB
lcube[iz,il,0]=phi_sm_2[il]*dv*0.125
;AGN same for now but fix 
lcube[iz,il,1]=phi_sm_2[il]*dv*0.125
endfor

endfor

end

pro plot_ccd

COMMON database,info,frame,lum_0,phi_0,dphi_0,zs,lcube,sedpar,lam,ext,powpah,filts,iter

wset, info.win_id5

data=LoadData(2)
white =GetColor('White',1)
black = GetColor('Black',2)
red = GetColor('Red',3)
blue = GetColor('Blue',4)

;plot,[0.0,20],[-1,35.0],color=white,/nodata,xtitle=TextoIDL('(F8/F4.5)'),ytitle=TeXtoIDL('(F24/F8.0)')

plot,[0.1,100],[0.1,15],color=white,/nodata,xtitle=TextoIDL('(F250/F24)'),ytitle=TeXtoIDL('(F8/F3.6)'),/xlog,/ylog,charsize=1.5,title='OUTPUT'
zs_old=zs

;input actual data for supersample
restore,'allfluxes_final.sav'
col1_data=f250/f24
col2_data=f8/f3
oplot,col1_data,col2_data,psym=3 ;,symthick=3

zs=zs_old

;plot,[-0.5,1.8],[-1.2,0.3],/nodata,xtitle=TextoIDL('log(F170/F70)'),ytitle=TeXtoIDL('log(F500/F250)')

;divide f8/f4.5 into 100 bins of 0.2 and f24/f8 in bins=0.5
colimsize=100
colim=fltarr(colimsize,colimsize)
binx=0.2
biny=0.5

area=2.5/1.d3 ;in sq.deg
stertodeg=(!pi/180.d0)^2.0

td=30.0
pfrac=0.6
tauv=10.0
tau2=25.0
alpha=4.0
alpha_hot=1.8
mstar=11.2

rms=[0,0,0.0032,0.0032,0,0.015,0.04,0,3.0,5.0,5.0,6.0,0] ;for the FLS [SoME numbers are arbitrary!, double check!]
zdet=[-1]

binf=4

col1_secondary=[-1]
col2_secondary=[-1]

for iz=1,n_elements(zs)-1 do begin ;omit z=0
    z=zs[iz]
    print,'z=',z
    for il=0,n_elements(lum_0)-5,4 do begin ;omit last 2 lbol as >13.5 and skip every second one for faster
        lbol=lum_0[il]
        num_sb=lcube[iz,il,0]*area*stertodeg ;absolute number of sources
        num_agn=lcube[iz,il,0]*area*stertodeg
        if(num_sb lt 1) then num_sb=num_sb*randomu(seed)
        if(num_sb lt 1) then goto,skipplot
        for in=long(0),num_sb,binf do begin
;            if(lbol lt 11.2) then mstar=lbol
;            if(lbol gt 11.2) then mstar=11.2
;            if(lbol lt 12.0) then tau2=5.0
;            if(lbol gt 12.4) then pfrac=0.3+0.2*randomn(seed)
;            if(lbol gt 12.0) then tau2=10.1+3.6*randomn(seed)

            f=sedmodel(z,lbol,mstar,tauv,alpha,tau2,pfrac,td)

;test whether or not detected and if so keep the values
;            if((f[3] gt 5.0*rms[3]) and (f[5] gt 5.0*rms[5]) and
;            (f[6] gt 5.0*rms[6])) then begin
;GO1 selection
;            if((f[6] ge 0.8) and (alog10((8.0/24.0)*f[6]/f[5]) ge 0.5) and (alog10((0.64/24.0)*f[6]/f[0]) ge 1.0)) then begin
;GO2 selection

;our specific 24um flux selection
            ;print,z,lbol,f[6]
             if(f[6] ge 0.8) then begin
                f[2]=f[2]+rms[2]*randomn(seed)
                f[3]=f[3]+rms[3]*randomn(seed)
                f[5]=f[5]+rms[5]*randomn(seed)
                f[6]=f[6]+rms[6]*randomn(seed)
;                f[8]=f[8]+rms[8]*randomn(seed)
                f[10]=f[10]+rms[10]*randomn(seed)
;                col1=(f[5]/f[3])
;                col2=(f[6]/f[5])
                col1=(f[10]/f[6])
                col2=(f[5]/f[3])
                if(iter eq 1) then oplot,[col1],[col2],psym=4,color=white
                if(iter eq 2) then oplot,[col1],[col2],psym=4,color=green
                if(iter eq 3) then oplot,[col1],[col2],psym=4,color=red

                col1_secondary=[col1_secondary,(f[5]/f[4])]
               col2_secondary=[col2_secondary,(f[6]/f[5])]

 ;               col2_secondary=[col2_secondary,(f[6]/f[0])]
                
                if(zdet[0] eq -1) then begin
                    zdet=z
                    cib_contr=f
                endif
                if(zdet[0] ne -1) then begin
                    zdet=[zdet,z]
                    cib_contr=cib_contr+f
                endif
                matchx=col1/binx
                matchy=(col2+10.0)/biny
                if((matchx ge 0) and (matchx lt colimsize) and (matchy ge 0) and (matchy lt colimsize)) then colim[matchx,matchy]=colim[matchx,matchy]+1
;                oplot,[col1],[col2],psym=3
;                print,z,lbol,num_sb,col1,col2
                
            endif else begin
                goto,skipplot
            endelse
        endfor
        skipplot:
    endfor
endfor

;help,colim
;tv,bytsc(colim)
area_st=area*stertodeg ;((!pi/180.d0)^2.0)

wset, info.win_id3
col1_data=f8/f4
col2_data=f24/(f8/1.d3)

;col2_data=f24/(fr/1.d3)

;print,fr
;print,col2_data
;print,col2_secondary

plot,col1_data,col2_data,yrange=[0,20],ystyle=1,psym=3,xtitle=TeXtoIDL('F8/F4.5'),ytitle=TeXtoIDL('F24/F8'),charsize=1.5,/nodata,color=white,title='OUTPUT'
oplot,col1_data,col2_data,psym=3,thick=5
oplot,col1_secondary,col2_secondary,psym=4,color=white

nzplot,zdet,binf
cibplot,(cib_contr/area_st),binf


end

PRO nzplot,zdet,binf
COMMON database, info, frame,lum_0,phi_0,dphi_0,zs,lcube,sedpar,lam,ext,powpah,filts,iter

data=LoadData(2)
white =GetColor('White',1)
black = GetColor('Black',2)
red = GetColor('Red',3)
blue = GetColor('Blue',4)

wset,info.win_id4

h=histogram(zdet,binsize=0.2,min=0,max=4.0,locations=xz)

ymax=max(h)
yr=[0.0,ymax*binf*1.4]
yr=[0.0,45]
if(iter eq 1) then plot,xz-0.1,h*binf,psym=10,xtitle=TeXtoIDL('redshift'),ytitle=TeXtoIDL('Number of sources'),color=white,linestyle=2,xrange=[0,4.0],yrange=yr,ystyle=1,xstyle=1,charsize=1.5,title='OUTPUT'

if(iter eq 2) then oplot,xz-0.1,h*binf,psym=10,color=green
if(iter eq 3) then oplot,xz-0.1,h*binf,psym=10,color=red

legend,position=[1.5,35],linestyle=[0,2],['observed','simulated'],color=white,box=0 ;,/nobox

;GO2 data
data_go2=read_ascii('supersample-withz',data_start=1)
zgo2=data_go2.field01[2,*]
hnz_go2=histogram(zgo2,binsize=0.2,min=0,max=4.0,locations=xz)
oplot,xz-0.1,hnz_go2,psym=10 ;,color=green

end

PRO cibplot,cib_contr,binf
COMMON database, info, frame,lum_0,phi_0,dphi_0,zs,lcube,sedpar,lam,ext,powpah,filts,iter

wset,info.win_id6

data=LoadData(2)
white =GetColor('White',1)
black = GetColor('Black',2)
red = GetColor('Red',3)
blue = GetColor('Blue',4)

bands=[0.64,0.79,3.6,4.5,8.0,24.0,70.0,100.0,170.0,250.0,500.0,850.0]
print,cib_contr[6] ;for testing purposes at 24um
cib_contr=cib_contr*3.d0/bands ;10^(-15)W/m2
cib_contr=cib_contr/(10.^6) ;to get it into cib units
cib_contr=cib_contr/20.0

plot,bands,cib_contr*binf,psym=10,xtitle=TeXtoIDL('\lambda [\mum]'),ytitle=TeXtoIDL('\nu I_{\nu} [nW m^{-2} st^{-1}]'),color=white,/xlog,/ylog,xrange=[0.1,1000.0],yrange=[0.1,5.d2],charsize=1.5,title='OUTPUT'

;Hauser & Dwek 2001 CIB data (only points that passed isotropy test
;and claim detection)
wav=[2.2,2.2,3.5,60,100,140,240,240,240]
pow=[23,20,12,28,25,32,17,14,13]

;oplot,wav,pow,psym=4,symsize=3,color=white

end


function sedmodel,agnfrac,lbol,mstar,tauv,alpha_hot,alpha,tau2,pfrac,td,lum
COMMON database, info, frame,lum_0,phi_0,dphi_0,zs,lcube,sedpar,lam,ext,powpah,filts,iter

z=2.0

conv = 2480.25*1.d0*4.d0*!pi
const=4799.0
freq=3.d0/lam

;fluxes=[0,0,0,0,0]

sspages = [0.0055,0.0255,0.1005,0.2905,0.6405,0.9005,1.405,2.505,5.05,11.05]

t=ageuniv(z)

good = where(sspages lt t[0],count)

age = sspages[good[n_elements(good)-1]]

if(age lt 1.0) then sspfile = 'ssp_'+gstring(fix(1.d3*age))+'Myr_z02.spec'

if((age ge 1.0) and (age lt 5.0)) then sspfile = 'ssp_'+gstring(fix(age))+'.'+gstring(10.*(age-fix(age)))+'Gyr_z02.spec'

if(age ge 5.0) then sspfile = 'ssp_'+gstring(fix(age))+'Gyr_z02.spec'

ssp = read_ascii(sspfile,data_start=6)

lam_ssp = ssp.field1[0,*]/1.d4
flux_ssp = ssp.field1[1,*]

pow_ssp = (flux_ssp*3.0/lam_ssp)
pow_ssp = pow_ssp*(10.^(mstar)/pow_ssp[6542]) ;scale to 1.6um [check!]
pow_ssp_new = interpolate(pow_ssp,findex(lam_ssp,lam))

powstars=pow_ssp_new

;warm component as a collection of temperatures [K]
;tcont=[1500,1000,800,600,400,300,200,150,100,80]
norm=1.0/alog(10.0)
beta=1.5
td=sedpar.td

aw=(1.0/1.d2) ;0.001*(10.^a[1])

!except=0 ;turn off exceptions
;warm = (aw*(freq^(1.d0-alpha))*exp(-tau2*ext))
;cold=(1.d11)*(freq^(4.0+1.5))/(exp(const*freq/td)-1)
;hf=where((lam gt 3) and (lam lt 25))

;at high pah fractions assume fixed ratio
;if(pfrac lt 0.5) then begin
;    lcold=total(cold)
;    lwarm=total(warm[hf])
;    sffrac=lcold/(lwarm+lcold)
;    sffrac_want=pfrac
;    aw=2.0*sffrac/sffrac_want
;    warm=aw*warm
;    lwarm=total(warm[hf])
;    sffrac=lcold/(lwarm+lcold)
;endif

;freq=3.0/lam

;parms=[?,?,ac,td,ah,aw,
    ;ac=0.01*10.^(parms[2])
;    ah=0.000001*10.^(parms[4])
;    aw=0.000001*10.^(parms[5])
    freqo=0.15 ;parms[10]
    fohot=0.08 ;parms[11]
    ;alpha=parms[7]
    beta=2.0

    cold=(freq^4.0)*(1.d0-exp(-((freq/freqo)^beta)))/(exp(const*freq/td)-1.d0)
    hot=freq/((freq/fohot)^(alpha)*exp(0.5*freq)+(freq/fohot)^(-0.5)+(freq/(0.3*fohot))^(-3.0))
    warm=(freq^(1.0-3.0))*exp(-0.27/freq)

;for id=0,n_elements(tcont)-1 do begin
;    ac=(10.^(2.9*alpha-3.34)*aw*tcont[id]^(1.0-alpha))*((const/tcont[id])^(4.d0+beta))/(10.d0^(0.69*(1.0+beta)))
;    xd=const*(freq)/tcont[id]
;    dust=ac*(freq^(4.0+beta)/(exp((xd) > (-1e-38) < 1e38)-1.0)) ;to avoid overflow messages
;    if(id eq 0) then dusttot=dust
;    if(id gt 0) then dusttot=dusttot +dust
;endfor

;warm=dusttot*exp(-tau2*ext)

d = distance(z)

;ycont=warm+cold
ycont=warm+1.d5*cold

lcont = (d^2.0)#(ycont*conv)
good=where((lam gt 3) and (lam lt 1000))
ltot=alog10(total(lcont[good])*0.005*alog(1.d1))
ac=10.^(lbol-ltot)


    ah=3.0*agnfrac*(total(ac*cold))/(total(hot*exp(-tau2*ext)))
    ycont=aw*warm+ah*hot*exp(-tau2*ext)+ac*5.d9*cold


good = where((lam ge 5) and (lam le 15))
ratio = total(ac*ycont[good])/total(powpah[good])
apah=ratio*5.0*pfrac/(1.-pfrac)
astars=1.0

;mfit = astars*powstars+apah*powpah+aw*warm+ah*hot*exp(-tau2*ext)+ac*5d10*cold
mfit=ah*hot*exp(-tau2*ext)+apah*powpah+aw*warm+ac*5.d9*cold

lum = (d^2.0)#(mfit*conv)

;lum = lum+pow_ssp_new*exp(-tauv*ext) ;(1.-exp(-tauv*ext_new))/tauv
lum=lum+0.05*powstars*exp(-2.0*ext)

		sh = fix(alog10(1.0+z)/0.005 + 0.5)
		d=distance(z)
		nufnu = 3.0*(1.0/d^2.0)#(lum/conv)
		nufnuobs = shift(nufnu,sh)
                ;print,'max fnuobs',max(nufnuobs)
fr = ((0.64/2.99)*total(nufnuobs*filts.r)/total(filts.r))
fi = ((0.79/2.99)*total(nufnuobs*filts.i)/total(filts.i))
f3 = ((3.6/2.99)*total(nufnuobs*filts.irac1)/total(filts.irac1))
f4 = ((4.5/2.99)*total(nufnuobs*filts.irac2)/total(filts.irac2))
f5 = ((5.8/2.99)*total(nufnuobs*filts.irac3)/total(filts.irac3))
f8 = ((8.0/2.99)*total(nufnuobs*filts.irac4)/total(filts.irac4))
f24 = ((24.0/2.99)*total(nufnuobs*filts.mips1)/total(filts.mips1))
f70 = (70.0/2.99)*total(nufnuobs*filts.mips2)/total(filts.mips2)
f100 = (100.0/2.99)*total(nufnuobs*filts.pacs100)/total(filts.pacs100)

f170 = ((170.0/2.99)*total(nufnuobs*filts.s170)/total(filts.s170))
f250 = (250.0/2.99)*total(nufnuobs*filts.spire250)/total(filts.spire250)
f350 = (350.0/2.99)*total(nufnuobs*filts.spire250)/total(filts.spire350)
f500 = (500.0/2.99)*total(nufnuobs*filts.spire250)/total(filts.spire500)
;f250 = (250.0/2.99)*nufnuobs[640]
;f500 = (500.0/2.99)*nufnuobs[700]

f850 = ((850.0/2.99)*total(nufnuobs*filts.s850)/total(filts.s850))

;imag = -2.5*alog10(fi)+16.57-0.45

fluxes=[fr,fi,f3,f4,f5,f8,f24,f70,f100,f170,f250,f350,f500,f850]

skipsed:

return,fluxes

end

PRO seds
COMMON database, info, frame,lum_0,phi_0,dphi_0,zs,lcube,sedpar,lam,ext,powpah,filts,iter

wset,info.win_id2

data=LoadData(2)
white =GetColor('White',1)
black = GetColor('Black',2)
red = GetColor('Red',3)
blue = GetColor('Blue',4)

lum=lam

fluxes=sedmodel(sedpar.agnfrac,sedpar.lbol,sedpar.mstar,sedpar.tauv,sedpar.alpha_hot,sedpar.alpha,sedpar.tau2,sedpar.pfrac,sedpar.td,lum)

if(iter eq 1) then plot,lam,lum,color=white,/xlog,/ylog,ystyle=1,xstyle=1,xrange=[0.3,950],yrange=[10.^9,7*10.^13],xtitle=TeXtoIDL('rest wavelength [\mum]'),ytitle=TeXtoIDL('Luminosity'),charsize=1.5,title='INPUT'

if(iter eq 2) then oplot,lam,lum,color=green
if(iter eq 3) then oplot,lam,lum,color=red

end
