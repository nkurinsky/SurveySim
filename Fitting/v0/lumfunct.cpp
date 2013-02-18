#include "lumfunc.h" //See this header for more information on the following funcitons

void lumfunct::set_phi0(double val){
  phi0=val;
}

void lumfunct::set_L0(double val){
  L0=val;
}

void lumfunct::set_alpha(double val){
  alpha=val;
}

void lumfunct::set_beta(double val){
  beta=val;
}

void lumfunct::set_p(double val){
  p=val;
}

void lumfunct::set_q(double val){
  q=val;
}

//this function returns the number of sources (per Mpc^3) for a given L-z pair
double lumfunct::get_nsrcs(double redshift,double lum){
  double nsrcs;
  double t1=pow(10,phi0)*pow((1.+redshift),p);
  double t2=pow(10,L0)*pow((1.+redshift),q);
  double ratio,lf;
  ratio=pow(10,lum)/t2;
  lf=t1/(pow(ratio,alpha)+pow(ratio,beta));
  //cout<<lf<<endl;
  //nsrcs=static_cast<int>(lf+0.5);
  nsrcs=lf;
  return nsrcs;
}

