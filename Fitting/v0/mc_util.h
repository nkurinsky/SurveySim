/*
  Utility Classes for Monte Carlo Fitter
  Written by Noah Kurinsky, 10/31/13

  These classes are used for controlling acceptance probabilities of chain steps
  and tracking accepted moves for m chains
*/

class MetropSampler {
 private:
  int nchains;
  long *accept_total;
  long *iteration_total;
  double tmax;
  bool accepted;
  double ideal_acceptance;
 public:
  MetropSampler(int nchains,double maxTemp, double idealpct);
  bool accept(int chainnum, double de, double normtime);
  double acceptance(int chainnum);
  double mean_acceptance();
  bool anneal();
  ~MetropSampler();
}

class MCChains {
 private:
  int allwidth;
  int nchains;
  int npars;
  int chainwidth;
  int nruns;
  int i;
  double chi_min;  
  double *bestpars;
  int *chainlength;
  valarray<valarray<double>> chains;
 public:
  MCChains(int nchains, int npars, int nruns);
  bool add_link(int chain, double pars, double chisqr);
  bool converged();
  bool save(string filename);
  ~MCChains();
}
