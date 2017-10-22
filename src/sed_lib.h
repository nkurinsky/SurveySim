// -*-c++-*-

#ifndef SED_LIB_H
#define SED_LIB_H

#include "functions.h"
#include "filters.h"
#include "cosmo.h"
#include "alg_lib/interpolation.h"
#include <map>
#include <tuple>

//for filter convolution integrals
#define SL_INT_SIZE 1000
#define SL_INT_PRECISION 1E-2
#define FILTER_NUM 3

using namespace CCfits;

class sed{
 private:
  gsl_interp_accel *acc;
  gsl_spline *spline;
  bool interp_init;
  alglib::spline2dinterpolant s;
  bool zdim;
  double brange[2];
  double zrange[2];
public:
  sed();
  sed(double *f,double *bands, int bnum);
  sed(double *f,double *bands, int bnum, double *z, int znum);
  double get_flux(double band, double redshift=0.0);
  ~sed();
};

class sed_lib{
 private:
  std::vector<unique_ptr<sed> > seds;
  unique_ptr<double[]> lums;
  unsigned int lnum;
  unsigned int tnum;
  bool redshiftq, typeq;
  //for interpolation
  gsl_interp_accel *acc;
  gsl_spline *spline;
  vector<alglib::spline2dinterpolant> flux_interpolator;
  bool interp_init;
  int interp_znum;
  int logflag;
  double interp_zmin;
  double interp_dz;
  //SEDs and evolution
  double brange[2];
  //filter convolution and interpolation (for speed)
  filter_lib filters;
  vector<string> types;
  void initialize_filter_fluxes(int logflag);
  double interpolate_flux(double lum, double redshift, short sedtype, short filter_id);
  double convolve_filter(short lum_id, double redshift, short sedtype, short filter_id);
  friend double flux_yield(double wavelength, void * params);
  gsl_integration_workspace *w;
 public:
  sed_lib(string fitsfile, int nz, double zmin, double dz);
  bool load_filters(string fitsfile,int logflag);
  void get_filter_info(string names[], double limits[], double errors[], double skewErrs[]);
  double get_flux(double lum, double redshift, short sedtype, double band);
  double get_filter_flux(double lum, double redshift, short sedtype, short filter_id);
  int get_lnum() {return lnum;}
  unordered_map<string, double> get_fracs();
  double get_dl();
  void get_lums(double luminosities[]);
  int get_tnum() {return tnum;}
  int get_snum() {return seds.size();}
  ~sed_lib();
};

double flux_yield(double wavelength, void * params);

struct flux_yield_params{
  sed *SED;
  filter *FILT;
  double z;
  flux_yield_params(sed *mySED, filter *myFILT, double redshift);
};


#endif
