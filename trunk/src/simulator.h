// -*-c++-*-

#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "obs_lib.h"
#include "sed_lib.h"
#include "lumfunct.h"
#include "hist_lib.h"
#include "functions.h"
#include "constants.h"
#include "cosmo.h"

//Storage structure for each individual source
struct sprop{
  //source properties (as determined by distributions)
  double redshift;
  double luminosity;
  double weight;
  double fluxes[3];
  double c1,c2;
  sprop();
  sprop(double z,double f[],double lum, double w); //constructor initializes variables
  //as well as automatically computes colors if it is a valid operation
  friend ostream &operator<<(ostream &out,sprop c);
};

struct products{
  double chisqr;
  valarray<double> dndz;
  valarray<double> dnds[3];
  products() : chisqr(0){
    dndz.resize(1);
  }
  products(int nz, int ns[]);
};

class simulator{
private:
  products last_output;
  vector<sprop> sources;
  //bool simulated;
  lumfunct *lf;
  std::unique_ptr<sed_lib> seds;
  std::unique_ptr<obs_lib> observations;
  std::unique_ptr<hist_lib> diagnostic;
  std::unique_ptr<numberCounts> counts[3];
  string filters[3];
  double band_errs[3];
  double flux_limits[3];
  double color_exp;
  double area;
  double dz;
  double zmin;
  int nz;
  int ns;
  string filterFile;
  RandomNumberGenerator rng;
  void initialize_filters();
  void initialize_counts();
 public:
  simulator();
  simulator(string filterfile, string obsfile, string sedfile);
  bool load_filter_lib(string file);
  bool load_filter(short filt_id, string name, double error, double flim);
  void set_size(double area,double dz,double zmin,int nz,int ns);
  void set_color_exp(double val, double zcut=1000);
  void set_lumfunct(lumfunct *lf);
  void set_sed_lib(string sedfile);
  void set_obs(string obsfile);
  void reset();
  products simulate();
  double model_chisq() { return last_output.chisqr; }
  bool save(string outfile);
  ~simulator();
};

#endif
