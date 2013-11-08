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

bool MetropSample::anneal(){
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


//chains implementation
MCChains::MCChains(int nchains, int npars, int nruns){
  chainwidth = npars+1;
  allwidth = chainwidth*nchains;
  this->npars = npars;
  this->nruns = nruns;
  this->nchains = nchains;
  chains.resize(allwidth);
  for(i=0;i<allwidth;i++)
    chains[i].resize(nruns);
  bestpars = new double[npars];
  chainlength = new int[nchains];
  for (i=0;i<nchains;i++)
    chainlength[i]=0;
  chi_min = 1000.0;
}

bool MCChains::add_link(int chain, double pars[], double chisqr){
  static int cbase;
  
  if (chain >= nchains){
    printf("MCChains::ERROR: Invalid chain number in add_link");
    return false;
  }

  if (chisqr < chi_min){
    for (i=0;i<npars;i++)
      bestpars[i]=pars[i];
    chi_min = chisqr;
  }
  
  if (chainlength[chain] >= nruns){
    printf("MCChains::ERROR: Invalid chain length in add_link");
    return false;
  }
  
  cbase=chain*chainwidth;
  for (i=0;i<npars;i++)
    chains[cbase+i][chainlength[chain]] = pars[i];
  chains[cbase+npars][chainlength[chain]] = chisqr;
  
  chainlength[chain]++;
  return true;
}

bool MCChains::converged(){
  return false;
}

bool MCChains::save(string filename, string parnames[]){

  using namespace CCfits;
  std::auto_ptr<FITS> pFits(0);
  
  try{
    pFits.reset(new FITS(filename,Write));
  }
  catch (CCfits::FITS::CantOpen){
    std::cerr << "Unknown Error Occurred, Can't save chain" << endl;
    return false;
  }
  
  std::vector<string> colnames(allwidth,"CHISQ");
  std::vector<string> colunits(allwidth,"-");
  std::vector<string> colform(allwidth,"e13.5");
  string hname("Chain");
  string cnum;
  
  printf("Output Columns\n");
  for(int j=0;j<nchains;j++){
    sprintf(cnum,"%d",j);
    for (i=0;i<npar;i++){
      colnames[j*chainwidth+i] = parnames[i]+cnum;
    }
    colnames[j*chainwidth+npar] += cnum;
  }
  
  Table *newTable = pFits->addTable(hname,nruns,colnames,colform,colunits,AsciiTbl);
  
  for(i=0;i<allwidth;i++)
    printf("\t%s\n",colnames[i]);
  newTable->column(colnames[i]).write(chain[i],1);
  return true;
}

MCChains::~MCChains(){
  
}
