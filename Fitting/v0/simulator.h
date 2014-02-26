#include "obs_lib.h"
#include "sed_lib.h"
#include "lumfunct.h"
#include "hist_lib.h"

#define LUMPARS 7

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
  products(){
    dndz.resize(1);
  }
  products(int nz, int ns[]);
};

struct clementsBins{
  //temporary structure created to compute clements'10 style number counts
  int bnum[3];
  //bin minima [mJy]
  valarray<double> b250;
  valarray<double> b350;
  valarray<double> b500;
  //bin widths [Jy]
  valarray<double> db250;
  valarray<double> db350;
  valarray<double> db500;
  //bin minima [Jy]^2.5
  valarray<double> S250;
  valarray<double> S350;
  valarray<double> S500;
  clementsBins();
};

class simulator{
 private:
  products *last_output;
  clementsBins dndsInfo;
  int binNum(int band, double flux);
  vector<sprop> sources;
  //bool simulated;
  lumfunct *lf;
  sed_lib *seds;
  obs_lib *observations;
  hist_lib * diagnostic;
  double bands[3];
  double band_errs[3];
  double flux_limits[3];
  double color_exp;
  double area;
  double dz;
  double zmin;
  int nz;
  int ns;
 public:
  simulator(){
    int bnumtemp[] = {0,0,0};
    last_output = new products(0,bnumtemp);
    last_output->chisqr=0;
    color_exp=0; //default to no color evolution
    area = pow((M_PI/180.0),2.0); //default to 1sq degree
    zmin = 0.1; //default to 0.1-6.0, 0.1 steps
    nz = 59;
    dz = 0.1;}
  simulator(double b[],double b_err[],double f_lims[],string obsfile,string sedfile);
  void set_bands(double b[],double b_err[],double f_lims[]);
  void set_size(double area,double dz,double zmin,int nz,int ns);
  void set_color_exp(double val);
  void set_lumfunct(lumfunct *lf);
  void set_sed_lib(string sedfile);
  void set_obs(string obsfile);
  void reset();
  products simulate();
  double model_chisq() { return last_output->chisqr; }
  bool save(string outfile);
  ~simulator();
};
