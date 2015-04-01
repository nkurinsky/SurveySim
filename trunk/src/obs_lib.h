// -*-c++-*-

#ifndef OBS_LIB_H
#define OBS_LIB_H

#include "functions.h"

using namespace CCfits;

class obs{
 public:
  double fluxes[3];
  double c1,c2;
  obs();
  obs(double f[], axis_type axes[]);
  double get_flux(int band);
  void get_colors(double &c1,double &c2);
  ~obs();
};

class obs_lib{
 private:
  vector<obs*> observations;
  double zp[3]; //zero points for magnitude-flux conversions
 public:
  obs_lib(string fitsfile, axis_type axes[], double flim[]);
  double get_flux(int i,int band);
  void get_colors(int i,double &c1,double &c2);
  void get_all_colors(double* &c1,double* &c2);
  int get_snum(){
    return observations.size();}
  ~obs_lib();
};

#endif
