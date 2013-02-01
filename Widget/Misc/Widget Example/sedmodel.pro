function sedmodel,z,lbol,mstar,tauv,alpha,tau2,pfrac,td

restore,'intervectors.save'

conv = 2250.d0*4.d0*!pi
const=4799.0

sspages = [0.0055,0.0255,0.1005,0.2905,0.6405,0.9005,1.405,2.505,5.05,11.05]

t=ageuniv(z)

good = where(sspages lt t[0],count)

age = sspages[good[n_elements(good)-1]]

if(age lt 1.0) then sspfile = '/Users/annie/bc03/templates/ssp_'+gstring(fix(1.d3*age))+'Myr_z02.spec'

if((age ge 1.0) and (age lt 5.0)) then sspfile = '/Users/annie/bc03/templates/ssp_'+gstring(fix(age))+'.'+gstring(10.*(age-fix(age)))+'Gyr_z02.spec'

if(age ge 5.0) then sspfile = '/Users/annie/bc03/templates/ssp_'+gstring(fix(age))+'Gyr_z02.spec'

ssp = read_ascii(sspfile,data_start=6)

lam = ssp.field1[0,*]/1.d4
flux = ssp.field1[1,*]

pow = (flux*3.0/lam)
pow = pow*(10.^(mstar)/pow[6542])
lam_ssp = lam
pow_ssp = pow

bandobs = [0.64,0.79,3.6,4.5,5.8,8.0,24.0,70,170.0,850.0] ;R and I
bandrest = bandobs/(1.+z)

inps = read_ascii('inputs2.dat')
lglam=inps.field01[0,*]
ext = inps.field01[1,*]
powpah = inps.field01[3,*]

freq = 2.99/(10.^(lglam))
lam=10.^(lglam)

pow_ssp_new = interpolate(pow_ssp,findex(lam_ssp,lam_new))

warm = ((5.0/1.d2)*(freq^(1.d0-alpha))*exp(-tau2*ext))
cold=(1.d11)*(freq^(4.0+1.5))/(exp(const*freq/td)-1)
hf=where((lam gt 3) and (lam lt 25))

;at high pah fractions assume fixed ratio
if(pfrac lt 0.5) then begin
    lcold=total(cold)
    lwarm=total(warm[hf])
    sffrac=lcold/(lwarm+lcold)
    print,'SF frac: ',sffrac
    sffrac_want=pfrac
    aw=2.0*sffrac/sffrac_want
    warm=aw*warm
    lwarm=total(warm[hf])
    sffrac=lcold/(lwarm+lcold)
    print,'SF new: ',sffrac
endif

d = distance(z)

ycont=cold
ycont[hf]=warm[hf]
;scale l12 so get Lbol right
lcont = (d^2.0)#(ycont*conv)
good=where((lam gt 3) and (lam lt 1000))
ltot=alog10(total(lcont[good])*0.005*alog(1.d1))
ac=10.^(lbol-ltot)

;scale to get PAH-fraction right
good = where((lam ge 5) and (lam le 15))
ratio = total(ac*ycont[good])/total(powpah[good])
apah=ratio*pfrac/(1.-pfrac)

mfit = apah*powpah+ac*ycont
lum = (d^2.0)#(mfit*conv)

lum_new = interpolate(lum,findex(lam,lam_new))

lum_new = lum_new+pow_ssp_new*(1.-exp(-tauv*ext_new))/tauv

		sh = fix(alog10(1.0+z)/0.005 + 0.5)
		d=distance(z)
		nufnu = (1.0/d^2.0)#(lum_new/conv)
		nufnuobs = shift(nufnu,sh)

fr = abs((bandobs[0]/2.99)*total(nufnuobs*filtr)/total(filtr))
fi = abs((bandobs[1]/2.99)*total(nufnuobs*filti)/total(filti))
f3 = abs((bandobs[2]/2.99)*total(nufnuobs*filt3)/total(filt3))
f4 = abs((bandobs[3]/2.99)*total(nufnuobs*filt4)/total(filt4))
f5 = abs((bandobs[4]/2.99)*total(nufnuobs*filt5)/total(filt5))
f8 = abs((bandobs[5]/2.99)*total(nufnuobs*filt8)/total(filt8))
f24 = abs((bandobs[6]/2.99)*total(nufnuobs*filt24)/total(filt24))

f70 = abs((bandobs[7]/2.99)*total(nufnuobs*filt70)/total(filt70))
f170 = abs((bandobs[8]/2.99)*total(nufnuobs*filt170)/total(filt170))
f850 = abs((bandobs[9]/2.99)*total(nufnuobs*filt850)/total(filt850))

f250 = (250.0/2.99)*nufnuobs[640]
f500 = (500.0/2.99)*nufnuobs[700]
;col1 = alog10((8.d0/24.d0)*f24/f8)
;col2 = alog10((bandobs[0]/bandobs[3])*f24/fr)

;imag = -2.5*alog10(fi)+16.57-0.45

;return,[f24,col1,col2,imag,f3]
return,[fr,fi,f3,f4,f5,f8,f24,f70,f170,f250,f500,f850]

end
