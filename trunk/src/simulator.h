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
  double flux_limits[3];
  double color_exp;
  double area;
  double dz;
  double zmin;
  int nz;
  string obsFile;
  string modelFile;
  axis_type axes[2];
  RandomNumberGenerator rng;
  void initialize_filters();
  void initialize_counts();
  long num_sources(double z, double l, double dl);
 public:
  simulator();
  simulator(const Configuration &config);
  simulator(string modelfile, string obsfile, string sedfile, axis_type axes[]);
  bool load_filters(string file);
  void set_diagnostic_xaxis(axis_type option);
  void set_diagnostic_yaxis(axis_type option);
  void set_diagnostic_axes(axis_type xopt, axis_type yopt);
  void set_size(double area,double dz,double zmin,int nz);
  void set_color_exp(double val, double zcut=1000);
  void set_lumfunct(lumfunct *lf);
  void set_sed_lib(string sedfile);
  void set_obs(string obsfile);
  void reset_obs();
  void reset();
  products simulate();
  double model_chisq() { return last_output.chisqr; }
  bool save(string outfile);
  ~simulator();  
};

#endif
