//Noah Kurinsky
//7/3/2012
/* This header file contains the prototypes for model and model_lib
classes. Model class is used to store and manipulate one model instance.
Model_lib class stores and manipulates all model classes for all models.
Model_lib also builds itself from a FITS file containing model and 
information about the parameters it contains
 */

#include "functions.h" //header containing standard includes

using namespace CCfits; //namespace specific to FITS reading

// class for storing an individual model, all public as it lies within
// the model_lib class and has no permissions but to itself
class model{
 private:
  gsl_interp_accel *acc;
  gsl_spline *spline;
 public:
  vector<double> fluxes;
  //vector<double> model_params;
  model(double *bands,double *fluxes,int fluxnum);
  //this function will redshift the desired flux and extrapolate a flux density
  double get_flux(double band,double redshift);
  ~model();
};

// class for managing model instances
class model_lib{
 private:
  model ** models;
  double *bands;
  int mod_num;
  int band_num;
  int pnum;
  //dynamic memory variable (for arbitrary pnum)
  int * psize;
  gsl_interp_accel * acc;
  gsl_interp *interp;
 public:
  //constructor reads in models from FITS file, initializes all variables
  //to match the properties of the models
  model_lib(string fitsfile); 
  //gets a redshifted flux from a given model number
  double get_flux(double ids[],double band,double redshift);
  //procedures for converting between model index, parameter specific id,
  //and parameter values corresponding to those ids
  int index(int ids[]);
  double index(double ids[]);
  int * get_ids(int index);
  //returns number of model parameters
  int get_pnum(){
    return pnum;}
  //returns number of values for each parameter (ps should have dim=pnum)
  void get_psize(int ps[]);
  ~model_lib();
};
