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

class sed_lib{
 private:
  vector<sed*> seds;
  vector<double> lums;
  double *bands;
  int bandnum;
 public:
  sed_lib(string fitsfile);
  double get_flux(double lum,double band);
  int get_snum(){
    return seds.size();}
  ~sed_lib();
};
