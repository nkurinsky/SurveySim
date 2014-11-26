pro update_wise

data=mrdfits('test/test.fits')
help,data,/struct
f1=data[0,*]
f2=data[1,*]
f3=data[2,*]
f4=data[3,*]
print,min(f1),min(f2),min(f3),min(f4)
plot,f1,f2,psym=3,/xlog,/ylog
header=headfits('test/test.fits')

fxbhmake,bhdr,176,'Data'
sxaddpar,bhdr,'NAXIS1',48
sxaddpar,bhdr,'TFIELDS',8 ;4 fluxes+4 error

sxaddpar,bhdr,'TTYPE1','F3.4'
sxaddpar,bhdr,'TTYPE2','F4.6'
sxaddpar,bhdr,'TTYPE3','F12'
sxaddpar,bhdr,'TTYPE4','F22'
sxaddpar,bhdr,'TTYPE5','E3.4'
sxaddpar,bhdr,'TTYPE6','E4.6'
sxaddpar,bhdr,'TTYPE7','E12'
sxaddpar,bhdr,'TTYPE8','E22'

sxaddpar,bhdr,'TUNIT1','mJy'
sxaddpar,bhdr,'TUNIT2','mJy'
sxaddpar,bhdr,'TUNIT3','mJy'
sxaddpar,bhdr,'TUNIT4','mJy'
sxaddpar,bhdr,'TUNIT5','mJy'
sxaddpar,bhdr,'TUNIT6','mJy'
sxaddpar,bhdr,'TUNIT7','mJy'
sxaddpar,bhdr,'TUNIT8','mJy'

sxaddpar,bhdr,'TFORM1','D'
sxaddpar,bhdr,'TFORM2','D'
sxaddpar,bhdr,'TFORM3','D'
sxaddpar,bhdr,'TFORM4','D'
sxaddpar,bhdr,'TFORM5','D'
sxaddpar,bhdr,'TFORM6','D'
sxaddpar,bhdr,'TFORM7','D'
sxaddpar,bhdr,'TFORM8','D'

sxaddpar,bhdr,'F1COL','F3.4'
sxaddpar,bhdr,'F2COL','F4.6'
sxaddpar,bhdr,'F3COL','F12'
sxaddpar,bhdr,'F4COL','F22'

sxaddpar,bhdr,'F1MIN',0.02,'Minimum flux1, these data'
sxaddpar,bhdr,'F2MIN',0.03,'Minimum flux2, these data'
sxaddpar,bhdr,'F3MIN',2.4,'Minimum flux3, these data'
sxaddpar,bhdr,'F4MIN',63.0,'Minimum flux4, these data'
sxaddpar,bhdr,'F1ERR',0.0067,'rms1, these data 3sig=min_f'
sxaddpar,bhdr,'F2ERR',0.02,'rms2, these data 3sig=min_f'
sxaddpar,bhdr,'F3ERR',0.8,'rms3, these data 3sig=min_f'
sxaddpar,bhdr,'F4ERR',21.,'rms4, these data 3sig=min_f'

;print,bhdr

;tmp=mrdfits('/Users/annie/students/noah_kurinsky/SurveySim/IDL/observation.fits',1,htmp)
;print,htmp

fxbcreate,5,'tmp.fits',bhdr
print,bunit ;,extension
;print,fxbcolnum(5,'F1COL')
fxbwritm,bunit,[1,2,3,4],data

fxbfinish,bunit

;fxwrite,'wise.fits',data,header
end
