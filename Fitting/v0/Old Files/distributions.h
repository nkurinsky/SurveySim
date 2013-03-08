/*
  Noah Kurinsky
  8/2/2012
  This header file provides the interface for the classes that manage both
  the redshift and model distributions
*/

#include "functions.h"

using namespace std;

//Redshift distribution class
class distribution{
 private:
  vector<double> values;
 public:
  //Constructors initialize either flat or gaussian distribution or fixed value
  distribution(gsl_rng * r,int size,double mean,double sigma,double range[]);
  distribution(gsl_rng * r,int size,double range[]);
  distribution(int size,double value);
  //Functions for adding to existing distribution
  void add(gsl_rng * r,int size,double mean,double sigma,double range[]);
  void add(gsl_rng * r,int size,double range[]);
  void add(int size,double value);
  //Functions which return samples or properties of the distribution
  int size(){
    return values.size();}
  double get(int i);
  void get_all(double retvals[]);
};

//Used by Lumdist as function/parameters passed to GSL Integrator
struct lf_params{double OmegaM;double Lambda0;};
double lf(double x,void * params);

//Luminosity function class
class lumfunct{
 private:
  bool init; //is the function initialized?
  double phi0; //Scale at z=0
  double L0;  //Location of the knee at z=0
  //Function Slope Parameters
  double alpha;
  double beta;
  double p,q; //Parameters which evolve Phi and L
  double exp_range[2]; //logarithmic range of L/L0 ratios
  double flux_limit;  //flux limit from sample for 350 microns
  double zmax;
  double H0,Lambda0,OmegaM; //Cosmological Parameters
  double ** luminosities; //dynamic 2D array for storing numerical luminosities
  //gsl variables for random luminosity generation and indexing
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_ran_discrete_t ** g_array;
  //luminosity distance calculator
  double lumdist(double z);
  gsl_interp_accel *ldist_acc;
  gsl_spline *ldist_spline;
  double *ldists;
  double *z_ldists;
 public:
  lumfunct(double flimit); //Sets default parameters, initializes function
  //For changing parameters
  void set_phi0(double val);
  void set_L0(double val);
  void set_alpha(double val);
  void set_beta(double val);
  void set_p(double val);
  void set_q(double val);
  void set_flux_limit(double val);
  void set_range(double lower,double upper);
  void set_H0(double val);
  void set_Lambda0(double val);
  void set_OmegaM(double val);
  //For extracting flux scale for given redshift
  double get_amp(double redshift, double &luminosity);
  //Initializes function, called from constructor or individually
  void initialize();
  ~lumfunct();
};

