#include "agn_frac.h"

//Fraction of galaxies as a function of luminosity and redshift that are AGN (note ALL AGN including obscured and unobscured, can try to sub-divide at some future point.
double agn_frac::get_agn_frac(double lum, double redshift){

  //the Chen,Hickox et al. 2013 relation between average(Log(Lx)) and Log(Lir)
  //here Lx is in ergs/s while lum (i.e. Lir) is in Lsun;
  float av_lglx=30.37+1.05*lum;
  
  //the Lutz et al. 2004 relation between Log(Lx) and Log(L6um) for Seybert 1s
  //both luminosities as in ergs/s;
  float av_lgl6um=av_lglx-0.41;

  //the average lgl6 to lg(Lagn_ir) conversion including converting to Lsun
  // for now a placeholder only, will need to calculate properly
  float av_lgAGNir=av_lgl6um-7-26-log10(3.86)+0.4;

  //the logic here is that if the average Lir_agn is comparable to the Lir of the given bin, essentially 100% of the galaxies are AGN and scales from there. With the Chen+13 relation, we end up with fairly low AGN fractions even at the highest luminosities.
  return pow(10,(av_lgAGNir-lum));
}
