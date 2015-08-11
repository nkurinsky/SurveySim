// -*-c++-*-

#ifndef AGN_FRAC_H
#define AGN_FRAC_H

#include <math.h>
#include "lumfunct.h"
#include "functions.h"

class agn_frac{
 private:
  lumfunct *lf;
  int _types;
  bool generate;
  bool hasComposites;
  RandomNumberGenerator rng;
  map <tuple<double,double>,double> values;
  map <double,double> lumPower;
  map <double,double> evolZ;
  double _t1;
  double _t2;
  double _fagn0;
  double _zbt;
 public:
  agn_frac(int agn_types);
  void set_lumfunct(lumfunct *lf);
  void set_params(double lpars[]);
  void set_t1(double t1);
  void set_t2(double t2);
  void set_fagn0(double fagn0);
  void set_zbt(double zbt);
  double get_agn_frac(double lum, double redshift);
  int get_sedtype(double lum, double redshift);
};


#endif
