#include "numberCounts.h"

NumberCounts::NumberCounts(const string &name) :  _name(name), _initialized(false), _verbose(false){

}

NumberCounts::NumberCounts(const valarray<double> &fluxes, const double area, const string &name) : NumberCounts(name) {
  initialize(fluxes,area);
}

bool NumberCounts::initialize(const valarray<double> &fluxes, const double area, const string &name){

  if(fluxes.min() > 0){
    double mean,N,sigma;
    double diff,nbins,binlow;
    valarray<double> logf(log10(fluxes));
    
    _range[0] = logf.min();
    _range[1] = logf.max();
    
    N=static_cast<double>(logf.size());
    mean = logf.sum()/N;
    sigma = sqrt((pow(logf-mean,2)).sum()/(N-1.0));
    
    _dS = 3.49*sigma/pow(N,0.33333333);
    nbins = ceil((_range[1]-_range[0])/_dS);
    _nbins = static_cast<int>(nbins);
    
    diff = _dS*nbins-(_range[1]-_range[0])/2.0;
    _range[0] = _range[0] - diff;
    _range[1] = _range[1] + diff;
    
    _bin_center.resize(_nbins,0.0);
    _scale_factors.resize(_nbins,1.0);
    for(int i=0.0;i<_nbins;i++){
      binlow = static_cast<double>(i)*_dS+_range[0];
      _bin_center[i] = binlow+_dS/2;
      _scale_factors[i] = pow(pow(10,_bin_center[i]),2.5)/(pow(10,binlow)*(pow(10,_dS)-1));
    }
    
    compute(fluxes,area,_counts);

    if(name != "")
      _name = name;

    _initialized = true;
  }
  else{
    _initialized = false;
    printf("NumberCounts::initialized Error: Invalid flux found in input array\n");
  }
  
  return _initialized;
}

const string& NumberCounts::name() const{
  return _name;
}

void NumberCounts::setName(string &name){
  _name = name;
}

void NumberCounts::compute(const valarray<double> &fluxes, const double area, valarray<double> &counts){
  static int _range_violations(0);
  
  if(fluxes.min() < _range[0] or fluxes.max() > _range[1]){
    _range_violations++;
    if(_verbose or _range_violations == 1){
      printf("NumberCounts::compute Warning: Input fluxes exceed counts histogram range\n");
      if(not _verbose)
	printf("NumberCounts::compute \tFurther Warnings will be Suppressed\n");
    }
  }
  
  counts.resize(_nbins,0.0);
  int j;
  for(unsigned int i = 0; i < fluxes.size(); i++){
    j = static_cast<int>( ceil( ( log10(fluxes[i]) - _range[0] ) / _dS ) );
    if(j >= 0 and j < counts.size())
      counts[j]++;
    else
      if(_verbose)
	printf("ERROR: Count flux out of range\n");
  }
  counts *= _scale_factors/area;
}

void NumberCounts::setVerbose(bool flag){
  _verbose=flag;
}


const valarray<double>& NumberCounts::counts() const{
    return _counts;
}

const valarray<double>& NumberCounts::bins() const{
  return _bin_center;
}
