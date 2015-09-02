//Noah Kurinsky
//5/25/2014

#include "cosmo.h"

double lumdist(double zmax) {

  double sum,z,dz,ez;
  static std::unordered_map<double,double> lumdists;
  if(lumdists.count(zmax) == 0){
    dz=0.001;
    sum=0.0;
    for (z=(dz/2.0);z<=zmax;z+=dz){
      ez=sqrt(OL+OM*pow((1.0+z),3));
      sum+=1/ez;
    } 
    lumdists[zmax] = (DH*(1.0+zmax)*dz)*sum;
  }

  return lumdists[zmax];
}

//assumes area never changes in a given code iteration
double dvdz (double zmax, double area)
{
  static std::unordered_map<double,double> dvdzStore;
  if(dvdzStore.count(zmax) == 0){
    double dz2 = pow((1+zmax),2);
    double dl=lumdist(zmax);
    double ez=sqrt(OL+(OM*dz2*(1+zmax)));
    dvdzStore[zmax]=(DH*dl*dl*area)/(ez*dz2);
  }

  return dvdzStore[zmax];
}

