//Noah Kurinsky
//7/10/2012
//This header file provides the interface for the class which manages the model
//and redshift distributions as well as simulates the result of these distributions

#include "distributions.h"
#include "model_lib.h"

//Storage structure for each individual source
struct sprop{
  //source properties (as determined by distributions)
  double mod_num;
  double redshift;
  double luminosity;
  double bands[3];
  double fluxes[3];
  double c1,c2;
  sprop();
  sprop(double m,double z,double b[],double f[],double lum); //constructor initializes variables
  //as well as automatically computes colors if it is a valid operation
  friend ostream &operator<<(ostream &out,sprop c);
};

//simulation class, stores a model and redshift distribution, ensures that they are always identical dimensions
//[may be expanded to include storage for various models read in from FITS file]
class sim{
 private:
  int pnum;
  distribution* z;
  distribution** mods;
  model_lib* lib;
  vector<sprop> sources;
  bool simulated;
  bool initialized;
  gsl_rng * r;
  const gsl_rng_type * T;
  void init_rand(); //Initializes GNU random number generator (see functions.h)
 public:
  sim(model_lib* mods);
  void initialize_dists(double mean[],double sigma[],int size,double min[],double max[]);
  void add_dists(double mean[],double sigma[],int size,double min[],double max[]);
  //functions which return samples or properties of the distributions
  void add_noise(double sigma[]);
  sprop get(int i);
  void get(int i,double &rz,double m[]);
  void get(int i,double &rz,double &m);
  void get_color(int i,double &c1,double &c2);
  void get_all_colors(double *&c1,double *&c2);
  int size();
  //will eventually run the simulation for the established distributions
  bool simulate(double bands[]);
  bool simulate(double bands[],double bfluxerr[],double flim[],lumfunct &lums);
  //as most variables are dynamic, an explicit deconstructor is necessary
  ~sim();
};
