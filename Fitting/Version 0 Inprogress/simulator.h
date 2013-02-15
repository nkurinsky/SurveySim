#include "sim.h"
#include "obs_lib.h"
#include "hist_lib.h"

class simulator{
 private:
  sim *distributions;
  lumfunct *luminosity_function;
  model_lib *models;
  obs_lib *observations;
  hist_lib * diagnostic;
  double bands[3];
  double band_errs[3];
  double flux_limits[3];
  double distribution_size;
  double redshift_range[2];
 public:
  simulator();
  simulator(double b[],double b_err[],double f_lims[]);
  simulator(lumfunct *lf,string modfile,string obsfile);
  simulator(double b[],double b_err[],double f_lims[],lumfunct *lf,string modfile,string obsfile);
  void set_bands(double b[],double b_err[],double f_lims[]);
  void set_lumfunct(lumfunct *lf);
  void set_model_lib(string modfile);
  void set_obs_lib(string obsfile);
  void set_zrange(double zrange[]);
  void set_simulated_size(double size);
  double simulate(const gsl_vector *v,void *params);
  bool save(string outfile);
  ~simulator();
};

struct mparam{
  int p_id;
  double mean;
  double sigma;
  double min;
  double max;
};

struct fixed_params{
  simulator * sim;
  int pnum;
  int znum;
  mparam *p;
};

double simulate(const gsl_vector *v,void *params);
