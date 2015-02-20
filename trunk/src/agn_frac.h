// -*-c++-*-

#ifndef AGN_FRAC_H
#define AGN_FRAC_H

#include <math.h>
#include "lumfunct.h"

class agn_frac{
 private:
  lumfunct *lf;
 public:
  void set_lumfunct(lumfunct *lf);
  double get_agn_frac(double lum, double redshift);
  double get_agn_frac2(double lum, double redshift, int agntype);
};


#endif
