#include "agn_frac.h"
#include "lumfunct.h"

agn_frac::agn_frac(int agn_types){
  if(agn_types <= 1){
    _types=1;
    generate=false;
    printf("No AGN Detected\n");
  }
  else{
    _types=2;
    generate=true;
  }
  lf=NULL;
}

void agn_frac::set_lumfunct(lumfunct *lf){
  if(lf != NULL){
    this->lf = lf;
  }
  else
    cout << "ERROR: NULL Pointer Passed to Simulator" << endl;
}
//Fraction of galaxies as a function of luminosity and redshift that are AGN (note ALL AGN including obscured and unobscured, can try to sub-divide at some future point.
double agn_frac::get_agn_frac(double lum, double redshift){
  tuple<double,double> point(lum,redshift);
  
  if(values.count(point) == 0){

    //the Chen,Hickox et al. 2013 relation between average(Log(Lx)) and Log(Lir)
    //here Lx is in ergs/s while lum (i.e. Lir) is in Lsun;
    float av_lglx=30.37+1.05*lum;
    
    //the Lutz et al. 2004 relation between Log(Lx) and Log(L6um) for Seybert 1s
    //both luminosities as in ergs/s;
    float av_lgl6um=av_lglx-0.41;
    
    //the average lgl6 to lg(Lagn_ir) conversion including converting to Lsun
    // for now a placeholder only, will need to calculate properly
    static float offset(-7-26-log10(3.86)+0.4);
    float av_lgAGNir=av_lgl6um+offset;

  //the logic here is that if the average Lir_agn is comparable to the Lir of the given bin, essentially 100% of the galaxies are AGN and scales from there. With the Chen+13 relation, we end up with fairly low AGN fractions even at the highest luminosities.

    //store to speed up later calculations
    values[point]= pow(10,(av_lgAGNir-lum));
  }
  
  return values[point];
}

//here agntype = 0 (all agn) or 1 or 2 (includes Type 2 and reddened Type 1)
double agn_frac::get_agn_frac2(double lum, double redshift, int agntype){

    //Using the quasar LFs from Lacy et al. 2015, note that lum,lstar0,lstar, and phistar are all given in log_10
    //note the luminosities in that paper are given in terms of the 5um monochromatic restframe luminosity [ergs*Hz]
    
    float phistar,lstar0,gamma1,gamma2,k1,k2,k3;
    
    if (agntype == 0){
      phistar=-4.75;
      gamma1=1.07;
      gamma2=2.48;
      lstar0=31.92;
      k1=1.05;
      k2=-4.71;
      k3=-0.034;
    }

    if (agntype == 1){
      phistar=-5.18;
      gamma1=0.25;
      gamma2=2.68;
      lstar0=31.99;
      k1=0.537;
      k2=-5.48;
      k3=0.768;
    }
    
    if (agntype == 2){
      phistar=-4.98;
      gamma1=1.09;
      gamma2=2.61;
      lstar0=31.91;
      k1=1.165;
      k2=-4.45;
      k3=-0.23;
    }
    
    float eps=log10((1+redshift)/(1+2.5));
    float lstar=lstar0+k1*eps+k2*pow(eps,2)+k3*pow(eps,3);
    
    //the lgl5 [ergs/Hz] to lgLagn_ir [Lsun] conversion
    float lum5um=lum+33+log10(3.86)-13.7782-0.48;
    
    float denom=pow(pow(10,lum5um)/pow(10,lstar),gamma1)+pow(pow(10,lum5um)/pow(10,lstar),gamma2);
    float phi=pow(10,phistar)/denom;
    
    float phi_all,fagn;
    phi_all = lf->get_nsrcs(redshift,lum);
    
    //the AGN fraction is the ratio of AGN Luminosity Function computed here and the overall luminosity function (from lumfunc.cpp), ensure that this does not exceed 100%

    fagn=pow(10,(phi-phi_all));
    if (fagn > 1)
      fagn=1.0;

    return fagn;
    
}

int agn_frac::get_sedtype(double lum, double redshift){
  if(generate)
    return get_agn_frac(lum,redshift) > rng.flat(0,1) ? 1 : 0;
  else
    return 0;
}
