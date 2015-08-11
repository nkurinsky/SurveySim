// -*-c++-*-

#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "agn_frac.h"
#include "obs_lib.h"
#include "sed_lib.h"
#include "lumfunct.h"
#include "hist_lib.h"
#include "functions.h"
#include "constants.h"
#include "cosmo.h"
#include "numberCounts.h"
#include "simulator_utils.h"

class simulator{
private:
  products last_output;
  vector<sprop> sources;
  //bool simulated;
  lumfunct *lf;
  std::unique_ptr<agn_frac> fagns;
  std::unique_ptr<sed_lib> seds;
  std::unique_ptr<obs_lib> observations;
  std::unique_ptr<hist_lib> diagnostic;
  std::unique_ptr<NumberCounts> counts[3];
  string filters[3];
  double band_errs[3];
  double skew_errs[3];
  bool hasSkewErr[3];
  double flux_limits[3];
  double color_exp;
  double area;
  double dz;
  double zmin;
  int nz;
  int logflag;
  string obsFile;
  string modelFile;
  axis_type axes[2];
  RandomNumberGenerator rng;
  long num_sources(double z, double l, double dl);
  void reset();
  bool simflag; //if true, observations ignored
  void initial_simulation();
 public:
  simulator(const Configuration &config);
  void configure(const Configuration &config);
  void set_color_exp(double val, double zcut=1000);
  void set_fagn_pars(double lpars[]);
  void set_lumfunct(lumfunct *lf);
  products simulate();
  double model_chisq() { return last_output.chisqr; }
  bool save(string outfile);
  ~simulator();  
};

#endif
