//created by Anna Sajina (09/27/12) 
//based on earlier code by Noah Kurinsky

#include "functions.h"

using namespace std;

//==============================================================================
//Luminosity function class
//******************************************************************************
class lumfunct {
 private:
  bool init; //is lumfunc initialized?
  double phi0; //Scale at z=0
  double L0;  //Location of the knee at z=0
  //Function Slope Parameters
  double alpha;
  double beta;
  double p,q; //Parameters which evolve Phi and L
 public:
  //lumfunct(double lpar[6]);
//For changing parameters
  void set_phi0(double val);
  void set_L0(double val);
  void set_alpha(double val);
  void set_beta(double val);
  void set_p(double val);
  void set_q(double val);
  //for getting the number of sources/Mpc^3 for a given L,z pair
  double get_nsrcs(double redshift,double lum);
  //~lumfunct();
};


