#include "numberCounts.h"

NumberCounts::NumberCounts(const string &name) :  _name(name), _initialized(false), _verbose(false), _range_violations(0){

}

NumberCounts::NumberCounts(const valarray<double> &fluxes, const double area, const string &name) : NumberCounts(name) {
  initialize(fluxes,area);
}

bool NumberCounts::initialize(const valarray<double> &fluxes, const double area, const string &name){

  if(fluxes.min() > 0){
    double mean,N,sigma;
    double diff,nbins,binlow;
    valarray<double> logf(log10(fluxes));
    
    _range[0] = logf.min()-0.3; //accept sources 0.5 times dimmer than minimum
    _range[1] = logf.max()+1.0; //accept sources 10 times brighter than maximum
    
    N=static_cast<double>(logf.size());
    mean = logf.sum()/N;
    sigma = sqrt((pow(logf-mean,2)).sum()/(N-1.0));
    
    _dS = 7.0*sigma/pow(N,0.33333333);
    nbins = ceil((_range[1]-_range[0])/_dS);
    _nbins = static_cast<int>(nbins);
    
    diff = (_dS*nbins-(_range[1]-_range[0]))/2.0;
    if(diff > 0){
      _range[0] = _range[0] - diff;
      _range[1] = _range[1] + diff;
    }    

    _bin_center.resize(_nbins,0.0);
    _scale_factors.resize(_nbins,1.0);
    for(int i=0.0;i<_nbins;i++){
      binlow = static_cast<double>(i)*_dS+_range[0];
      _bin_center[i] = binlow+_dS/2;
      _scale_factors[i] = pow(pow(10,_bin_center[i])/1e3,2.5)/((pow(10,binlow)*(pow(10,_dS)-1))/1e3);
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

void NumberCounts::compute(const valarray<double> &fluxes_nolog, const double area, valarray<double> &counts){
  valarray<double> fluxes(log10(fluxes_nolog));
  
  if(fluxes.min() < _range[0] or fluxes.max() > _range[1]){
    _range_violations++;
    if(_verbose and _range_violations == 1){
      printf("NumberCounts::compute: Input fluxes exceed counts histogram range\n");
      printf("NumberCounts::compute:      Further Warnings will be Suppressed\n");
      printf("NumberCounts::compute:      Range is set by observations, most likely this is due\n");
      printf("NumberCounts::compute:      to bright low-redshift sources or sparse observations\n");
    }
  }
  
  counts.resize(_nbins,EMPTYBIN);
  int j;
  for(unsigned int i = 0; i < fluxes.size(); i++){
    j = static_cast<int>( ceil((fluxes[i] - _range[0]) / _dS ) );
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
