// -*-c++-*-

#ifndef AGN_FRAC_H
#define AGN_FRAC_H

#include <math.h>

class agn_frac{
 public:
  double get_agn_frac(double lum, double redshift);
  double get_agn_frac2(double lum, double redshift, int agntype);
};


#endif
