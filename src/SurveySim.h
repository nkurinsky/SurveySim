#ifndef SURVEYSIM_H
#define SURVEYSIM_H

//#include "functions.h"
#include "simulator.h"
#include "mc_util.h"

class SurveySim {
 private:
  Configuration q;
  int logflag;
  RandomNumberGenerator rng;
  
  //general variable declarations
  bool accept;
  bool saved = false;
  unsigned long i, m;
  unsigned long pi;
  double trial;
  
  //this array holds LF params
  double lpars[LUMPARS];
  double lmin[LUMPARS];
  double lmax[LUMPARS];
  double parfixed[LUMPARS];
  
  unique_ptr<string[]> parnames;
  lumfunct lf;
  simulator survey;
  
  products output;
  //ResultChain final_counts;
  
  MCChains mcchain;
  MCChains burnchain;
  
  vector<double> means;
  double chi_min; 
  double tchi_min;
  ParameterSettings pset;
  vector<double> results;
  MetropSampler metrop;
  
  vector<vector<double>> ptemp;
  vector<vector<double>> pcurrent;
  double **setvars;
  
  void initializeChains();
  void burnIn();
  void fit();
  void runExtra();
  void bestFit();
  void runOneStep(int m);
  void cleanUp();
 public:
  SurveySim(int argc, char** argv);
  void run();
  void save();
};

#endif
