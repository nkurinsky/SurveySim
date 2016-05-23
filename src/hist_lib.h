/*
  Noah Kurinsky
  7/9/2012
  This header provides the interface for the class that creates and performs operations on
  the two dimensional histograms used to fit the model to the data.

  Interface to be expanded to compute statistics on histograms
*/

#ifndef HIST_LIB_H
#define HIST_LIB_H

#include "functions.h"

//Needed for array manipulation by the GSL statistics procedures
#define H_STRIDE 1

class hist_lib{
 private:
  //variables for storage of dynamically created histograms
  double ** model_hist;
  double ** obs_hist;
  double ** comparison_hist;
  //properties of the histograms, set internally
  double xrange[2];
  double yrange[2];
  double xbinsize;
  double ybinsize;
  double chisq;
  int xysize;
  int osize,msize;
  double m_exc;
 public:
  //constructors
  hist_lib(); //All set to null
  hist_lib(double obs_c1[],double obs_c2[],int obs_size); //only observation histogram created
  //Functions for passing data to class if not already done so by constructor
  void init_obs(double c1[],double c2[],int size);
  bool init_model(double c1[],double c2[],int size);
  //Functions for outputting properties of histograms
  double get_chisq();
  bool write_fits(string filename);
  //Deconstructor for dynamic histogram arrays
  ~hist_lib();
 private:
  //functions which produce histograms
  double ** get_hist(double c1[],double c2[],int cnum); //compute histogram, fitting binsize
  double ** compute_hist(double c1[],double c2[],int cnum); //compute histogram for given inputs
  double fit_err();
  double poiss_err(int x);
};

#endif
