/*
  Utility Classes for Monte Carlo Fitter
  Written by Noah Kurinsky, 10/31/13

  These classes are used for controlling acceptance probabilities of chain steps
  and tracking accepted moves for m chains
*/

#include "functions.h"
#include <deque>

#define RECENT_NUM 100

class MetropSampler {
 private:
  int nchains;
  long *accept_total;
  long *iteration_total;
  double *previous;
  double temp;
  bool accepted;
  double ideal_acceptance;
  deque<bool> recent;
  unsigned long recent_num;
  gsl_rng *rgen;
 public:
  MetropSampler(int nchains,double maxTemp, double idealpct, gsl_rng *rgen);
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
  int i;
  double Rmax;
  double alpha;
  double chi_min;  
  double *bestpars;
  int *chainlength;
  valarray<valarray<double> > chains;
 public:
  MCChains(int nchains, int npars, int nruns);
  bool add_link(int chain, double pars[], double chisqr);
  void get_best_link(double pars[], double chisqr);
  bool set_constraints(double Rmax, double alpha);
  bool converged();
  bool save(string filename, string parnames[]);
  ~MCChains();
};
