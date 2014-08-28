// -*-c++-*-

#ifndef SED_LIB_H
#define SED_LIB_H

#include "functions.h"
#include "filters.h"
#include "cosmo.h"
#include <map>
#include <tuple>

//for filter convolution integrals
#define SL_INT_SIZE 1000
#define SL_INT_PRECISION 1E-3

using namespace CCfits;

class sed{
 private:
  gsl_interp_accel *acc;
  gsl_spline *spline;
  bool interp_init;
  double *fluxes;
 public:
  sed();
  sed(double *f,double *bands, int bnum);
  double get_flux(double band);
  ~sed();
};

//eventually, multiple SEDs per luminosity will be used in 
//specified porportions. Currently, one SED per luminosity
class sed_lib{
 private:
  std::vector<sed*> seds;
  double *lums;
  unsigned int lnum;
  unsigned int bandnum;
  gsl_interp_accel *acc;
  gsl_spline *spline;
  double brange[2];
  double color_exp;
  filter_lib filters;
  double convolve_filter(short lum_id, double redshift, short filter_id);
  friend double flux_yield(double wavelength, void * params);
  gsl_integration_workspace *w;
 public:
  sed_lib(string fitsfile);
  bool init_filter_lib(string file);
  bool load_filter(short filter_id, string name);
  double get_flux(double lum, double redshift, double band);
  double get_filter_flux(double lum, double redshift, short filter_id);
  void set_color_evolution(double exp);
  int get_lnum(){
    return lnum;}
  double get_dl();
  void get_lums(double luminosities[]);
  int get_snum(){
    return seds.size();}
  ~sed_lib();
};

double flux_yield(double wavelength, void * params);

struct flux_yield_params{
  sed_lib * wlib;
  short lum_ind;
  short filt_ind;
  double z;
  flux_yield_params();
};


#endif