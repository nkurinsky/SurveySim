#include "functions.h"

using namespace CCfits;

class sed{
 private:
  gsl_interp_accel *acc;
  gsl_spline *spline;
  bool interp_init;
 public:
  double *fluxes;
  sed();
  sed(double *f,double *bands, int bnum);
  double get_flux(double band);
  ~sed();
};

//eventually, multiple SEDs per luminosity will be used in 
//specified porportions. Currently, one SED per luminosity
class sed_lib{
 private:
  vector<sed*> seds;
  double *lums;
  int lnum;
  int bandnum;
  gsl_interp_accel *acc;
  gsl_spline *spline;
 public:
  sed_lib(string fitsfile);
  double get_flux(double lum,double band);
  int get_lnum(){
    return lnum;}
  void get_lums(double luminosities[]);
  int get_snum(){
    return seds.size();}
  ~sed_lib();
};
