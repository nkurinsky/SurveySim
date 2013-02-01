FUNCTION ageuniv,z

OM=0.27
OL=0.73
Ho=71.
c=2.99e8
Omega=1.

bbage = 13.7 ;wmap

pi=!Pi

dz=0.001
dimen=N_ELEMENTS(z)
lback=findgen(dimen)*0.0

j=0
while (j lt dimen) do begin
maxz=(z[j]*1000.)

FOR i=1,maxz DO BEGIN
	zr=i/1000.
	lback[j]=lback[j]+dz/((1.+zr)*sqrt((1.+zr)*(1.+zr)*(OM*zr+1.)-zr*(2.+zr)*OL))

ENDFOR	

;need to get lback in 'Gyr' so conversion constant

conv=(1.d3/Ho)

lback[j]=lback[j]*conv

j=j+1

endwhile

RETURN,bbage-lback

END
