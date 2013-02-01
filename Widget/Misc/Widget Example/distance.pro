FUNCTION distance,z

;WMAP 7-year Komatsu et al. 2011
OM=0.272
OL=0.728
Ho=70.4

c=2.99e8
Omega=1.
pi=!Pi

dz=0.00001
dimen=N_ELEMENTS(z)
distance=findgen(dimen)*0.0

j=0
while (j lt dimen) do begin
maxz=(z[j]*100000.)

FOR i=long(1),maxz DO BEGIN
	zr=i/100000.
	distance[j]=distance[j]+dz/sqrt(OM*(1.+zr)^3.+OL)
ENDFOR	

;need to get distance in 'Mpc' so conversion constant

conv=c*(1.+z[j])/(Ho*1000.)

distance[j]=distance[j]*conv

j=j+1

endwhile

RETURN,distance

END
