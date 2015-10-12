#include "agn_frac.h"
#include "lumfunct.h"

agn_frac::agn_frac(int agn_types){
  hasComposites=false;
  hasCold=false;
  if(agn_types <= 1){
    _types=1;
    generate=false;
    printf("No AGN Detected\n");
  }
  else{
    _types=agn_types;
    generate=true;
    if (agn_types > 2)
      hasComposites=true;
    
    if(agn_types > 3)
      hasCold=true;

    if(agn_types > 4){
      printf("Warning: AGN Types larger than expected limit of 3\n");
      printf("Number of AGN Types: %d \n",agn_types);
    }
    
  }
  lf=NULL;

  _fagn0=0;
  _compFrac=0.0;
  _coldFrac=0.0;
  _t1=0;
  _t2=0;
  _zbt=0;
}

void agn_frac::set_lumfunct(lumfunct *lf){
  if(lf != NULL){
    this->lf = lf;
  }
  else
    cout << "ERROR: NULL Pointer Passed to Simulator" << endl;
}

void agn_frac::set_params(double lpars[]){
  set_t1(lpars[LF::parameter::t1]);
  set_t2(lpars[LF::parameter::t2]);
  set_zbt(lpars[LF::parameter::zbt]);
  set_fagn0(lpars[LF::parameter::fa0]);
  set_fComp(lpars[LF::parameter::fcomp]);
  set_fCold(lpars[LF::parameter::fcold]);
}

void agn_frac::set_t1(double t1){
  if(t1 != _t1){
    _t1=t1;
    evolZ.clear();
  }
}

void agn_frac::set_t2(double t2){
  if(t2 != _t2){
    _t2=t2;
    evolZ.clear();
  }
}

void agn_frac::set_fagn0(double fagn0){
  _fagn0=fagn0;
}

void agn_frac::set_zbt(double zbt){
  if(zbt != _zbt){
    _zbt=zbt;
    evolZ.clear();
  }
}

void agn_frac::set_fComp(double fComp){
  if(fComp < 0.0)
    _compFrac=0.0;
  else if(fComp > 1.0)
    _compFrac=1.0;
  else
    _compFrac=fComp;
}

void agn_frac::set_fCold(double fCold){
  if(fCold < 0.0)
    _coldFrac=0.0;
  else if(fCold > 1.0)
    _coldFrac=1.0;
  else
    _coldFrac=fCold;
}

double agn_frac::get_agn_frac(double lum, double redshift){

  //here fagn0 gives the fraction of AGN at L_ir=10^12 at z=0
  //zbt gives the break redshift for the sed type
  // t1 and t2 give the power relations below and after the break respectively

  if(evolZ.count(redshift) == 0){
    if(redshift <= _zbt) 
      evolZ[redshift]=pow((1+redshift),_t1);
    if(redshift > _zbt)
      evolZ[redshift]=pow((1+(redshift-_zbt)),_t1)*pow((1+redshift),_t2);
  }

  if(lumPower.count(lum) == 0){
    lumPower[lum]=pow((lum/12.000),12.000);
  }

  //the agn fraction for a given L,z bin
  double fagn=_fagn0*lumPower[lum]*evolZ[redshift];
  return (fagn > 1 ? 1.0 : fagn);
}

//assume half of the AGN are composites
//assume half of the SFG are cold 
int agn_frac::get_sedtype(double lum, double redshift){
  if(generate){
    float fagn, rnd;    
    fagn=get_agn_frac(lum,redshift);
    rnd=rng.flat(0,1);
    if(rnd > fagn){
      if(hasCold){
	rnd=rng.flat(0,1);
	return rnd > _coldFrac ? 0 : 3; 
      }
      else
	return 0;
    }
    else{
      if(hasComposites){
	rnd=rng.flat(0,1);
	return rnd > _compFrac ? 1 : 2;
      }
      else
	return 1;
    }
  }
  // return sfg 
  return 0;
}
