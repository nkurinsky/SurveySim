#include "obs_lib.h"
#include "sed_lib.h"
#include "lumfunct.h"
#include "hist_lib.h"

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
  vector<double> dndz;
  vector<double> dnds;
  products();
  products(int nz, int ns);
};

class simulator{
 private:
  double chisq;
  vector<sprop> sources;
  bool simulated;
  lumfunct *lf;
  sed_lib *seds;
  obs_lib *observations;
  hist_lib * diagnostic;
  double bands[3];
  double band_errs[3];
  double flux_limits[3];
  double distribution_size;
  void init_rand(); //Initializes GNU random number generator (see functions.h)
 public:
  simulator() {chisq=0;}
  simulator(double b[],double b_err[],double f_lims[],string obsfile,string sedfile);
  void set_bands(double b[],double b_err[],double f_lims[]);
  void set_lumfunct(lumfunct *lf);
  void set_sed_lib(string sedfile);
  void set_obs(string obsfile);
  void reset();
  products simulate(double area, int nz, double dz, double zmin, int ns);
  double model_chisq() { return chisq; }
  bool save(string outfile);
  ~simulator();
};
