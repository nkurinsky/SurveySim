/*
  Noah Kurinsky
  7/9/2012
  This header provides the interface for the class that creates and performs operations on
  the two dimensional histograms used to fit the model to the data.

  Interface to be expanded to compute statistics on histograms
*/

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
  double range[2];
  double binsize;
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
  bool init_model(double c1[],double c2[],double weights[],int size);
  //Functions for outputting properties of histograms
  double get_chisq();
  bool write_fits(string filename);
  //Deconstructor for dynamic histogram arrays
  ~hist_lib();
 private:
  //functions which produce histograms
  double ** get_hist(double c1[],double c2[],int cnum,double range[]); //compute histogram, fitting binsize
  double ** compute_hist(double c1[],double c2[],double weights[],int cnum); //compute histogram for given inputs
  double fit_err();
  double poiss_err(int x);
  double poiss[101] = {1.841,2.3,2.638,2.918,3.163,3.383,3.584,3.77,3.945,4.11,4.267,4.417,4.56,4.698,4.83,4.959,5.083,5.204,5.321,5.435,5.547,5.655,5.761,5.865,5.967,6.067,6.164,6.26,6.355,6.447,6.538,6.628,6.716,6.803,6.888,6.973,7.056,7.138,7.219,7.299,7.377,7.455,7.532,7.608,7.684,7.758,7.832,7.904,7.976,8.048,8.118,8.188,8.258,8.326,8.394,8.461,8.528,8.594,8.66,8.725,8.789,8.853,8.917,8.979,9.042,9.104,9.165,9.226,9.287,9.347,9.407,9.466,9.525,9.583,9.641,9.699,9.756,9.813,9.87,9.926,9.982,10.037,10.092,10.147,10.202,10.256,10.31,10.363,10.417,10.47,10.522,10.575,10.627,10.678,10.73,10.781,10.832,10.883,10.933,10.984,11.034};
};
