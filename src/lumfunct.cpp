#include "lumfunct.h" //See this header for more information on the following functions

lumfunct::lumfunct(LF::distribution dist){
  switch(dist){
  case LF::distribution::Schecter:
  case LF::distribution::DoublePowerLaw:
  case LF::distribution::ModifiedSchecter:
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
  case LF::parameter::p2:
    return p2;
  case LF::parameter::q2:
    return q2;
  case LF::parameter::zbp:
    return zbp;
  case LF::parameter::zbq:
    return zbq;
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
  case LF::parameter::p2:
    p2 = value;
    return;
  case LF::parameter::q2:
    q2 = value;
    return;
  case LF::parameter::zbp:
    zbp = value;
    return;
  case LF::parameter::zbq:
    zbq = value;
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
  lpars[LF::parameter::p2]    = p2;
  lpars[LF::parameter::q2]    = q2;
  lpars[LF::parameter::zbp]  = zbp;
  lpars[LF::parameter::zbq]  = zbq;
}

void lumfunct::set_params(double lpars[]){
  phi0  = lpars[LF::parameter::PHI0];
  L0    = lpars[LF::parameter::L0];
  alpha = lpars[LF::parameter::alpha];
  beta  = lpars[LF::parameter::beta];
  p     = lpars[LF::parameter::p];
  q     = lpars[LF::parameter::q];
  p2    = lpars[LF::parameter::p2];
  q2    = lpars[LF::parameter::q2];
  // break redshift in "p"
  zbp = lpars[LF::parameter::zbp];
  // break redshift in "q"
  zbq = lpars[LF::parameter::zbq];
}

//this function returns the number of sources (per Mpc^3 per dlogL) for a given L-z pair
double lumfunct::get_phi(double redshift,double lum){
  static double norm(pow(10,phi0));
  static double knee(pow(10,L0));
  static double lastphi0(phi0);
  static double lastL0(L0);

  static double zbptop(pow(1+zbp,p));
  static double zbqtoq(pow(1+zbq,q));
  static double lastp(p);
  static double lastzbp(zbp);
  static double lastq(q);
  static double lastzbq(zbq);

  if(phi0 != lastphi0){
    lastphi0=phi0;
    norm=pow(10,phi0);
  }
  if(L0 != lastL0){
    lastL0=L0;
    knee=pow(10,L0);
  }

  double t1(norm);
  double t2(knee);
  double ratio;
  
  if(redshift <= zbp) {
    t1*=pow((1.+redshift),p);
  }
  else{
    if((zbp != lastzbp) or (p != lastp)){
      lastzbp=zbp;
      lastp=p;
      zbptop = pow(1+zbp,p);
    }
    t1*=zbptop*pow((1.+redshift-zbp),p2);
  }

  if(redshift <= zbq) {
    t2*=pow((1.+redshift),q);
  } 
  else{
    if((zbq != lastzbq) or (q != lastq)){
      lastq=q;
      lastzbq=zbq;
      zbqtoq = pow(1+zbq,q);
    }
    t2*=zbqtoq*pow((1.+redshift-zbq),q2);
  }
  
  static unordered_map<int,double> luminosity;
  int mapped_lum = static_cast<int>(floor(lum*MULTIPLIER));
  if(luminosity.count(mapped_lum) == 0)
    luminosity[mapped_lum]=pow(10,lum);
  
  ratio=luminosity[mapped_lum]/t2;
  
  switch(_dist){
  case LF::distribution::Schecter:
    return t1*pow(ratio,alpha)*exp(-ratio);
  case LF::distribution::DoublePowerLaw:
    return t1/(pow(ratio,alpha)+pow(ratio,beta));
  case LF::distribution::ModifiedSchecter:
    return t1*pow(ratio,1-alpha)*exp(-1*pow(log10(1+ratio),2)/(2*pow(beta,2)));
  default:
    static int logflag=3;
    LOG_CRITICAL(printf("ERROR: Invalid LF dist (%i) requested\n",_dist));
    //constructor ensures we never get here, but it is necessary for compilation
    return 0;
  }
}
