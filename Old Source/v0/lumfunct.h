//created by Anna Sajina (09/27/12) 
//based on earlier code by Noah Kurinsky

#ifndef LUMFUNCT_H
#define LUMFUNCT_H

#include "functions.h"

using namespace std;

//==============================================================================
//Luminosity function class
//******************************************************************************
class lumfunct {
 private:
  //bool init; //is lumfunc initialized?
  double phi0; //Scale at z=0
  double L0;  //Location of the knee at z=0
  //Function Slope Parameters
  double alpha;
  double beta;
  double p,q; //Parameters which evolve Phi and L
  double zmax;
 public:
  //lumfunct(double lpar[7]); uncomment if you want a constructor
//For changing parameters
  void set_phi0(double val){phi0=val;}
  void set_L0(double val){L0=val;}
  void set_alpha(double val){alpha=val;}
  void set_beta(double val){beta=val;}
  void set_p(double val){p=val;}
  void set_q(double val){q=val;}
  void set_zmax(double val){zmax=val;}
  void set_params(double lpars[]);
  void get_params(double lpars[]);
  //for getting the number of sources/Mpc^3 for a given L,z pair
  double get_nsrcs(double redshift,double lum);
  //~lumfunct(); uncomment if you want a destructor
};

#endif
