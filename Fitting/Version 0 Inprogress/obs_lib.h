#include "functions.h"

using namespace CCfits;

class obs{
 private:
  gsl_interp_accel *acc;
  gsl_spline *spline;
  bool interp_init;
 public:
  double fluxes[3];
  double c1,c2;
  obs();
  obs(double *f,double *bands);
  double get_flux(double band);
  void get_colors(double &c1,double &c2);
  ~obs();
};

class obs_lib{
 private:
  vector<obs*> observations;
  double bands[3];
  double ferr[3];
  double flim[3];
 public:
  obs_lib(string fitsfile);
  double get_flux(int i,double band);
  void get_colors(int i,double &c1,double &c2);
  void get_all_colors(double* &c1,double* &c2);
  void info(double *b,double *flims,double *ferrs);
  int get_snum(){
    return observations.size();}
  ~obs_lib();
};
