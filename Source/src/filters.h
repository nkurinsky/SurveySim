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
  filter(string filtername, vector<double> band, vector<double> transmission);
  bool load(string filtername, vector<double> band, vector<double> transmission);
  double transmission(double wavelength);
  double low();
  double high();
  void print(bool all=false);
  ~filter();
 private:
  string name;
  double *lambda;
  double *response;x
  gsl_interp_accel *acc;
  gsl_spline *spline;
  bool init;
  double filter_limits[2];
  unsigned long filter_size;
  double trap_integrate(vector<double> lambda, vector<double> response);
};

class filter_lib{
 public:
  filter_lib();
  filter_lib(string ffilename);
  bool init(){
    return initialized;}
  bool load_library(string ffilename);
  bool load_filter(short num, string fname);
  filter& get(short num);
  void print();
 private:
  vector<filter_info> library;
  filter filters[3];
  filter dummy;
  bool initialized;
};

#endif
