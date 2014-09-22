//Noah Kurinsky -*-c++-*-
//7/5/2012
//This header file provides the interface for the various small functions
//used for constraining models

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

//standard includes used by all programs in simulation
#include <cstdio>
#include <stdio.h>
#include <string>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <memory>
#include <cctype>

//For FITS input/output
#include <CCFits/CCfits>

//GNU scientific library includes
// - for random distributions
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
// - for statistics procedures
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_statistics_int.h>
// - for interpolation procedures
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
// - for numerical integration
#include <gsl/gsl_integration.h>
// - for scientific constants
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>

#define LUMPARS 7
#define SDEG_PER_STER 3282.8

using namespace std;

enum axis_type{
  ColorF1F3,
  ColorF2F3,
  ColorF1F2,
  Flux1,
  Flux2,
  Flux3
};

axis_type set_axis_type(string &keyvalue);

double metric_value(const double& f1,const double &f2,const double &f3,const axis_type &opt);

string toLower(const string &oldstr);

class Configuration{
public:
  Configuration(int argc, char *argv[]);
  void print();
  double areaSteradian() const;
  void LFparameters(double lpars[],short type = 0);
  
  string outfile;
  string obsfile;
  string modfile;
  string sedfile;
  string filterfile;
  
  bool vary_cexp;
  bool oprint;

  int nz;
  int ns;
  int cind;
  unsigned long runs;
  unsigned long nchain;
  unsigned long burn_step;
  unsigned long burn_ratio;
  unsigned long conv_step;
  unsigned long burn_num;
  unsigned long nparams;

  double area;
  double dz;
  double zmax;
  double zmin;
  double rmax;
  double a_ci;
  double tmax;
  double tscale;
  double idealpct;
  double annrng;
  
  vector<int> param_inds;
  double LFParameters[LUMPARS][4];
  double colorEvolution[4];

  const short value = 0;
  const short fixed = 1;
  const short min = 2;
  const short max = 3;

  axis_type axes[2];

private:
  void load();
  
};

class RandomNumberGenerator{
public:
  RandomNumberGenerator();
  ~RandomNumberGenerator();
  double gaussian(double mean, double sigma, double min, double max);
  double flat(double min, double max);
private:
  const gsl_rng_type *T;
  gsl_rng *r;
};

#endif
