FUNCTION dvol,z,dz,domega

;domega is in steradians

	c=299792458.
	om=0.27
        ol=0.73
        ho=71. ;km/s/Mpc

        dl=distance(z)

        conv=c/(ho*1000.0)
;        a=((1.0+z)^2.0)*(om*z+1)-z*(2.0+z)*ol
        a=om*((1.0+z)^3.0)+ol
        dvol = conv*(dl^2.0)/(((1.0+z)^2.0)*sqrt(a))

        dvol=dvol*domega*dz ;in Mpc^3

;        dvol = dvol*(!pi/180.0)^2 ; in sq. deg.

RETURN,dvol

END
