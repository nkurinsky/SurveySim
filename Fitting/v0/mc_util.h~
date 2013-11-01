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
  int nchains;
  int npars;
  int nruns;
  double chi_min;
  
 private:
  
 public:
}
