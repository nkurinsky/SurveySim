// -*-c++-*-

#ifndef NUMBERCOUNTS_H
#define NUMBERCOUNTS_H

#include "functions.h"
#include <valarray>

using namespace std;

class NumberCounts {
  
public:
  NumberCounts(const string &name="NULL");
  NumberCounts(const valarray<double> &fluxes,
	       const double area,
	       const string &name="NULL");
  bool initialize(const valarray<double> &fluxes, 
		  const double area, 
		  const string &name="");
  void compute(const valarray<double> &fluxes, 
	       const double area, 
	       valarray<double> &counts);
  void setName(string &name);
  void setVerbose(bool flag=true);
  const string& name() const;
  const valarray<double>& counts() const;
  const valarray<double>& bins() const;
private:  
  valarray<double> _bin_center;
  valarray<double> _scale_factors;
  valarray<double> _counts;
  double _range[2];
  string _name;
  double _dS;
  int _nbins;
  bool _initialized;
  bool _verbose;
  int _range_violations;
};

#endif
