#include "lumfunct.h" //See this header for more information on the following funcitons

void lumfunct::get_params(double lpars[]){
  lpars[0] = phi0;
  lpars[1] = L0;
  lpars[2] = alpha;
  lpars[3] = beta;
  lpars[4] = p;
  lpars[5] = q;
  lpars[6] = zmax;
}

void lumfunct::set_params(double lpars[]){
  phi0 = lpars[0];
  L0 = lpars[1];
  alpha = lpars[2];
  beta = lpars[3];
  p = lpars[4];
  q = lpars[5];
  zmax = lpars[6];
}

//this function returns the number of sources (per Mpc^3) for a given L-z pair
//this is the main function which would needs to be changed in order to change the adopted luminosity function
//note need to also give it dlogl as the number of sources within a bin 
//changes depending on the width of the bin!
double lumfunct::get_nsrcs(double redshift,double lum){
  double nsrcs;
  double t1,t2;
  double ratio;

  string lf_form ("dbl_pl");
  //string lf_form='modSchecter';

  if(redshift <= zmax) {
    t1=pow(10,phi0)*pow((1.+redshift),p);
    t2=pow(10,L0)*pow((1.+redshift),q);
  } 
  if(redshift > zmax) {
    t1=pow(10,phi0)*pow((1.+zmax),p);
    t2=pow(10,L0)*pow((1.+zmax),q);
  }
  ratio=pow(10,lum)/t2;

  if (lf_form.compare("dbl_pl") == 0)
    nsrcs=t1/(pow(ratio,alpha)+pow(ratio,beta));
  if (lf_form.compare("modSchecter") == 0)
    nsrcs=t1*(pow(ratio,1-alpha)*exp(-(1./2*pow(beta,2.))*(log(1-ratio)**2.);

  return nsrcs;
}
