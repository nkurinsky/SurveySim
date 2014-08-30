//Anna Sajina
//10/16/2012
//This header keeps various cosmology functions

#ifndef COSMO_H
#define COSMO_H

#include "constants.h"
#include <math.h>
#include <unordered_map>

// luminosity distances
double lumdist(double zmax);
double dvdz (double zmax, double area);

#endif
