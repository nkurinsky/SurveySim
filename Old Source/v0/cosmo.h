//Anna Sajina
//10/16/2012
//This header keeps various cosmology functions

#ifndef COSMO_H
#define COSMO_H

#include "constants.h"

// luminosity distances
double lumdist(double zmax) {

  int i;
  double tmp,z1,d,dz,ez;

  dz=zmax/1000;
  tmp=0;
  z1=0;
  for (i=1;i<=1000;i++){
    z1+=dz;
    ez=sqrt(OL+OM*(1+z1)*(1+z1)*(1+z1));
    tmp+=dz/ez;
  } 

  d=d_hub*tmp*(1.+zmax);
  return (d);
}

double dvdz (double zmax, double area)
{
  static double ez;
  static double my_dvdz;
  static double dl;
  static double dz2;

  dz2 = pow((1+zmax),2);
  dl=lumdist(zmax);
  ez=sqrt(OL+(OM*dz2*(1+zmax)));
  my_dvdz=(d_hub*dl*dl*area)/(ez*dz2);

  return my_dvdz;
}

#endif
