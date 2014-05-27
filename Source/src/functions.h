//Noah Kurinsky
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

using namespace std;

//Generate random numbers by use of procedures outlined in GNU Scientific Library Ch. 18 (Sections 1-5)
double * random(gsl_rng * r,double range[],int size);

//For Generation of Gaussian random number generator distributions (GNU Scientific Library Ch 20.2)
double * gauss_random(gsl_rng * r,double range[],double mean,double sigma,int size);

//Generates spectral color of a galaxy from fluxes and bands
double get_color(double f1,double f2);

#endif
