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
  string filename;
};

class filter{
 public:
  filter();
  filter(string filtername, string filename);
  bool load(string filtername, string filename);
  double transmission(double wavelength);
  double low();
  double high();
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
};

class filter_lib{
 public:
  filter_lib();
  filter_lib(string ffilename);
  bool init(){
    return init;}
  bool load_library(string ffilename);
  bool load_filter(short num, string fname);
  filter& get(short num);
  void print();
 private:
  vector<filter_info> library;
  filter filters[3];
  filter dummy;
  bool init;
};

#endif
