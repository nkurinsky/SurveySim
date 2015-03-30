// -*-c++-*-

#ifndef SIMULATORUTILS_H
#define SIMULATORUTILS_H

#include "functions.h"

//Storage structure for each individual source
struct sprop{
  //source properties (as determined by distributions)
  double redshift;
  double luminosity;
  double fluxes[3];
  double c1,c2;
  sprop();
  //constructor initializes variables
  //as well as automatically computes colors if it is a valid operation
  sprop(double z,double f[],double lum, axis_type opts[]); 
  friend ostream &operator<<(ostream &out,sprop c);
};

struct products{
  double chisqr;
  valarray<double> dndz;
  valarray<double> dnds[3];
  products() : chisqr(0){
    dndz.resize(1);
  }
  products(int nz, int ns[]);
};

#endif
