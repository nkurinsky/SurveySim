//-*-c++-*-
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
  double p,q,p2,q2; //Parameters which evolve Phi and L
  double zbp,zbq;
  LF::distribution _dist;
 public:
  lumfunct(LF::distribution dist=LF::DoublePowerLaw);
  double get_param(LF::parameter p);
  void set_param(LF::parameter p, double value);
  void set_params(double lpars[]);
  void get_params(double lpars[]);
  double get_phi(double redshift,double lum);
};

#endif
