/*
  Comment
*/

#ifndef FILTERS_H
#define FILTERS_H

#include "functions.h"

#include <vector>
#include <cstdio>

struct filter_info{
  string name;
  string info;
  double scale;
  vector<double> band;
  vector<double> transmission;
};

class filter{
 public:
  filter();
  filter(string filtername, vector<double> band, vector<double> transmission,int logflag);
  bool load(string filtername, vector<double> band, vector<double> transmission,int logflag);
  double transmission(double wavelength);
  double low();
  double high();
  string get_name(){
    return name;
  }
  void print(bool all=false);
  ~filter();
 private:
  string name;
  double *lambda;
  double *response;
  gsl_interp_accel *acc;
  gsl_spline *spline;
  bool init;
  double filter_limits[2];
  unsigned long filter_size;
  double trap_integrate(vector<double> lambda, vector<double> response);
};

class filter_lib{
 public:
  int logflag;
  filter_lib();
  filter_lib(string fitsfile,int logflag);
  bool init(){
    return initialized;}
  bool load_filters(string fitsfile,int logflag);
  void filter_info(string names[], double fluxLimits[], double fluxErrors[], double skewErrors[]);
  filter& get(short num);
 private:
  filter filters[3];
  double limits[3];
  double errors[3];
  double skew_errors[3];
  filter dummy;
  bool initialized;
};

#endif
