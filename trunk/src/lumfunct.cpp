#include "lumfunct.h" //See this header for more information on the following funcitons

lumfunct::lumfunct(LF::distribution dist){
  switch(dist){
  case LF::distribution::DoublePowerLaw:
  case LF::distribution::ModifiedSchechter:
    _dist=dist;
    break;
  default:
    fprintf(stderr,"Invalid Distribution passed to Lumfunct constructor\n");
    _dist=LF::DoublePowerLaw;
  }
}

double lumfunct::get_param(LF::parameter par){
  switch(par){
  case LF::parameter::PHI0:
    return phi0;
  case LF::parameter::L0:
    return L0;
  case LF::parameter::alpha:
    return alpha;
  case LF::parameter::beta:
    return beta;
  case LF::parameter::p:
    return p;
  case LF::parameter::q:
    return q;
  case LF::parameter::zmax:
    return zmax;
  default:
    fprintf(stderr,"Error in lumfunct::get_param: Invalid parameter %i\n",par);
    return -999;
  }
}

void lumfunct::set_param(LF::parameter par, double value){
  switch(par){
  case LF::parameter::PHI0:
    phi0 = value;
    return;
  case LF::parameter::L0:
    L0 = value;
    return;
  case LF::parameter::alpha:
    alpha = value;
    return;
  case LF::parameter::beta:
    beta = value;
    return;
  case LF::parameter::p:
    p = value;
    return;
  case LF::parameter::q:
    q = value;
    return;
  case LF::parameter::zmax:
    zmax = value;
    return;
  default:
    fprintf(stderr,"Error in lumfunct::set_param: Invalid parameter %i\n",par);
    return;
  }
}

void lumfunct::get_params(double lpars[]){
  lpars[LF::parameter::PHI0]  = phi0;
  lpars[LF::parameter::L0]    = L0;
  lpars[LF::parameter::alpha] = alpha;
  lpars[LF::parameter::beta]  = beta;
  lpars[LF::parameter::p]     = p;
  lpars[LF::parameter::q]     = q;
  lpars[LF::parameter::zmax]  = zmax;
}

void lumfunct::set_params(double lpars[]){
  phi0  = lpars[LF::parameter::PHI0];
  L0    = lpars[LF::parameter::L0];
  alpha = lpars[LF::parameter::alpha];
  beta  = lpars[LF::parameter::beta];
  p     = lpars[LF::parameter::p];
  q     = lpars[LF::parameter::q];
  zmax  = lpars[LF::parameter::zmax];
}

//this function returns the number of sources (per Mpc^3) for a given L-z pair
//this is the main function which would needs to be changed in order to change the adopted luminosity function
//note need to also give it dlogl as the number of sources within a bin 
//changes depending on the width of the bin!
double lumfunct::get_nsrcs(double redshift,double lum){
  double nsrcs;
  double t1,t2;
  double ratio;
  
  if(redshift <= zmax) {
    t1=pow(10,phi0)*pow((1.+redshift),p);
    t2=pow(10,L0)*pow((1.+redshift),q);
  } 
  if(redshift > zmax) {
    t1=pow(10,phi0)*pow((1.+zmax),p);
    t2=pow(10,L0)*pow((1.+zmax),q);
  }
  ratio=pow(10,lum)/t2;
  
  switch(_dist){
  case LF::distribution::DoublePowerLaw:
    return t1/(pow(ratio,alpha)+pow(ratio,beta));
  case LF::distribution::ModifiedSchechter:
    return t1*pow(ratio,1-alpha)*exp(-1*pow(log(1-ratio),2)/(2*pow(beta,2)));
  default:
    //constructor ensures we never get here, but it is necessary for compilation
    return nsrcs;
  }
}
