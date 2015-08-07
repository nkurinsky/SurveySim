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

double agn_frac::get_agn_frac(double lum, double redshift, double fagn0,double zbt,double t1,double t2){

  //here fagn0 gives the fraction of AGN at L_ir=10^12 at z=0
  //zbt gives the break redshift for the sed type
  // t1 and t2 give the power relations below and after the break respectively
 
  double fagn; //the agn fraction for a given L,z bin
  double evol; //evolution term

  if(redshift <= zbt) 
    evol=pow((1+redshift),t1);
  if(redshift > zbt)
    evol=pow((1+(redshift-zbt)),t1)*pow((1+redshift),t2);

  fagn=fagn0*pow((lum/12.0),4)*evol;

  //ensure does not exceed maximum
    if (fagn > 1)
      fagn=1.0;

    return fagn;
}

//old version with AGN only
//int agn_frac::get_sedtype(double lum, double redshift){
//  if(generate)
//    return get_agn_frac(lum,redshift) > rng.flat(0,1) ? 1 : 0;
//  else
//    return 0;
//}

//new version, assume half of the AGN are composites
    int agn_frac::get_sedtype(double lum, double redshift,double fagn0,double zbt,double t1,double t2){
  float fagn,fcomp,fsfg,rnd;
  int stype;
  fagn=get_agn_frac(lum,redshift,fagn0,zbt,t1,t2);
  fcomp=fagn/2.0;
  fagn=fagn/2.0;
  fsfg=1-fagn-fcomp;
  //if(generate)
  rnd=rng.flat(0,1);
  stype= rnd < fsfg ? 0 : 1;
  if(stype == 1)
    stype=rnd < fsfg+fcomp ? 1 : 2;
    //e.g. fsfg=0.4, fcomp=0.3 and fagn=0.3
    //rnd<0.4 -> stype=0 i.e. sfg
    //rnd < 0.7 (0.4+0.3) stype=1 i.e. composite
    //rnd > 0.7 stype = 2 i.e. agn
  return stype;
}
