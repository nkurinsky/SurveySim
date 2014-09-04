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
  bool getCounts(valarray<double> &counts) const;
  bool getBins(valarray<double> &bins) const;
  void getName(string &name) const;
  void setName(string &name);
  void setVerbose(bool flag=true);
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
};

#endif
