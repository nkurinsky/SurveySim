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
  RandomNumberGenerator rng;
  map <tuple<double,double>,double> values;
 public:
  agn_frac(int agn_types);
  void set_lumfunct(lumfunct *lf);
  double get_agn_frac(double lum, double redshift,double fagn0,double zbt,double t1,double t2);
  int get_sedtype(double lum, double redshift,double fagn0,double zbt,double t1,double t2);
};


#endif
