/*
  Utility Classes for Monte Carlo Fitter
  Written by Noah Kurinsky, 10/31/13

  These classes are used for controlling acceptance probabilities of chain steps
  and tracking accepted moves for m chains
*/

#ifndef MC_UTIL_H
#define MC_UTIL_H

#include "functions.h"
#include <deque>

#define RECENT_NUM 100

class ParameterSettings{
 public:
  ParameterSettings(size_t nparams);
  void set(short pnum, double Minimum, double Maximum, double standardDeviation, double bestValue);
  vector<double> min;
  vector<double> max;
  vector<double> sigma;
  vector<double> best;
};

class MetropSampler {
 private:
  int nchains;
  long *accept_total;
  long *iteration_total;
  double *previous;
  double temp;
  bool accepted;
  double ideal_acceptance;
  double accept_buffer;
  double accept_ratio;
  deque<bool> recent;
  deque<double> trials;
  unsigned long recent_num;
  RandomNumberGenerator rng;
 public:
  MetropSampler(int nchains,double maxTemp, double idealpct, double acpt_buf);
  bool accept(int chainnum, double trial);
  double acceptance(int chainnum);
  double mean_acceptance();
  double acceptance_rate();
  bool anneal();
  void reset();
  ~MetropSampler();
};

class MCChains {
 private:
  int allwidth;
  int nchains;
  int npars;
  int chainwidth;
  int nruns;
  int convruns;
  int i;
  double Rmax;
  double alpha;
  double chi_min;  
  double *bestpars;
  int *chainlength;
  valarray<valarray<double> > chains;
  valarray<valarray<double> > rvals;
  valarray<valarray<double> > accepted;
 public:
  MCChains(int nchains, int npars, int nruns, int convstep);
  bool add_link(int chain, double pars[], double chisqr, bool accpt);
  void get_best_link(double pars[], double chisqr);
  void get_fit_results(double pars[], double sigma[]);
  void get_stdev(double sigma[]);
  bool set_constraints(double Rmax, double alpha);
  bool converged();
  bool save(string filename, string parnames[]);
  ~MCChains();
};

class ResultChain {
 private:
  unsigned long i;
  vector<vector<valarray<double> > > results;
  vector<double> chisqrs;
 public:
  ResultChain(int num_arrays, int nresults);
  bool add_link(valarray<double> arrays[], double chisqr);
  bool save(string filename, string resnames[]);
};

#endif
