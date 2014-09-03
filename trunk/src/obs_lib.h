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
  obs(double f[]);
  double get_flux(int band);
  void get_colors(double &c1,double &c2);
  ~obs();
};

class obs_lib{
 private:
  vector<obs*> observations;
  string filter[3];
  double ferr[3];
  double flim[3];
 public:
  obs_lib(string fitsfile);
  double get_flux(int i,int band);
  void get_colors(int i,double &c1,double &c2);
  void get_all_colors(double* &c1,double* &c2);
  void info(string filters[],double flims[],double ferrs[]);
  int get_snum(){
    return observations.size();}
  ~obs_lib();
};

#endif
