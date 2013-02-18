//#include "obs_lib.h"
#include "model_lib.h"
//#include "sed_models.h"
#include "lumfunc.h"
#include "hist_lib.h"

//Storage structure for each individual source
struct sprop{
  //source properties (as determined by distributions)
  int mod_num;
  double redshift;
  double luminosity;
  double bands[3];
  double fluxes[3];
  double c1,c2;
  sprop();
  sprop(int m,double z,double b[],double f[],double lum); //constructor initializes variables
  //as well as automatically computes colors if it is a valid operation
  friend ostream &operator<<(ostream &out,sprop c);
};

class simulator{
 private:
  double m_chi2;
  //sim *distributions;
  vector<sprop> sources;
  //bool simulated;
  lumfunct *luminosity_function;
  model_lib *models;
  //  obs_lib *observations;
  hist_lib * diagnostic;
  double bands[3];
  double band_errs[3];
  double flux_limits[3];
  double distribution_size;
  double redshift_range[2];
  void init_rand(); //Initializes GNU random number generator (see functions.h)
 public:
  simulator()
    {
      m_chi2=0;
    }; //this defines the constructor simulator
  simulator(double b[],double b_err[],double area,double f_lims[],double lpars[],string modfile,string obsfile,string sedfile);
  //sprop get(int i);
  void set_bands(double b[],double b_err[],double f_lims[]);
  void set_lumfunct(lumfunct *lf);
  void set_model_lib(string modfile);
  void set_obs_lib(string obsfile);
  void set_zrange(double zrange[]);
  void set_simulated_size(double size);
  void reset();
  double simulate(const gsl_vector *v,void *params);
  double model_chi2() { return m_chi2; }
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
