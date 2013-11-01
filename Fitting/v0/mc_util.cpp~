/*
  Utility Classes Function Definitions for Monte Carlo Fitter
  Written by Noah Kurinsky, 10/31/13

  These classes are used for controlling acceptance probabilities of chain steps
  and tracking accepted moves for m chains
*/

#include "mc_util.h"

MetropSampler::MetropSampler(int nchains, double maxTemp, double idealpct){
  this->nchains = nchains;
  accept_total = new long[nchains];
  iteration_total = new long[nchains];
  for (int i=0;i<nchains;i++){
    accept_total[i] = 0;
    iteration_total[i] = 0;
  }
  tmax = maxTemp;
  ideal_acceptance = idealpct;
}

bool MetropSampler::accept(int chainnum, double de, double normtime){
  static double itemp,tester,tmc;
  tmc = normtime*tmax;
  tester=de/tmc;
  
  //to avoid errors from taking e^x where x is very large
  if(tester <= -2){
    accepted = (de<0.0) ? true : false;
  }
  else{
    itemp=gsl_rng_uniform(r);   //want uniform random number from 0-1
    accepted = ((de < 0.0) or (itemp < exp(-de/tmc))) ? true : false;
  }
  
  iteration_total[chainnum]++;
  if (accepted)
    accept_total[chainnum]++;
  return accepted;
}

double MetropSample::acceptance(int chainnum){
  return (double(accept_total[chainnum])/(double(iteration_total[chainnum])));
}

double MetropSample::mean_acceptance(){
  double acc_temp;
  for (int i=0;i<nchains;i++){
    acc_temp += acceptance(i);
  }
  return acc_temp/double(nchains);
}

bool MetropSample::anneal(double idealpct){
  if(mean_acceptance() > ideal_acceptance){
    temp--;
    return false;
  }
  return true;
}

MetropSampler::~MetropSampler(){
  delete accept_total;
  delete iteration_total;
}
