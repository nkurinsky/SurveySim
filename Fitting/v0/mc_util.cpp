/*
  Utility Classes Function Definitions for Monte Carlo Fitter
  Written by Noah Kurinsky, 10/31/13

  These classes are used for controlling acceptance probabilities of chain steps
  and tracking accepted moves for m chains
*/

#include "mc_util.h"

int dcomp(const void * a,const void * b){
  
  if((*(double*)a) > (*(double*)b)) return 1;
  if((*(double*)a) < (*(double*)b)) return -1;
  else return 0;
}

MetropSampler::MetropSampler(int nchains, double maxTemp, double idealpct, gsl_rng *rgen){
  this->nchains = nchains;
  accept_total = new long[nchains];
  iteration_total = new long[nchains];
  previous = new double[nchains];
  for (int i=0;i<nchains;i++){
    accept_total[i] = 0;
    iteration_total[i] = 0;
    previous[i] = 1.0E+4;
  }
  temp = maxTemp;
  ideal_acceptance = idealpct;
  this->rgen = rgen;
}

bool MetropSampler::accept(int chainnum, double trial){
  static double itemp,tester;
  tester=(trial-previous[chainnum])/temp;
  
  //to avoid errors from taking e^x where x is very large
  if(tester <= -2){
    accepted = (tester<0.0) ? true : false;
  }
  else{
    itemp=gsl_rng_uniform(rgen);   //want uniform random number from 0-1
    accepted = ((tester < 0.0) or (itemp < exp(-tester))) ? true : false;
  }
  
  iteration_total[chainnum]++;
  if (accepted){
    accept_total[chainnum]++;
    previous[chainnum] = trial;
  }
  
  return accepted;
}

double MetropSampler::acceptance(int chainnum){
  return (double(accept_total[chainnum])/(double(iteration_total[chainnum])));
}

double MetropSampler::mean_acceptance(){
  double acc_temp = 0;
  for (int i=0;i<nchains;i++){
    acc_temp += acceptance(i);
  }
  return acc_temp/double(nchains);
}

bool MetropSampler::anneal(){
  if(mean_acceptance() > ideal_acceptance){
    temp -= 5;
    return true;
  }
  return false;
}

MetropSampler::~MetropSampler(){
  delete[] accept_total;
  delete[] iteration_total;
  delete[] previous;
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
  //defaults to 95% CI
  Rmax = 1.05;
  alpha = 0.05;
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

bool MCChains::set_constraints(double Rmax, double alpha){

  if (Rmax <= 1.0 or Rmax > 2.0){
    printf("MCChains::ERROR: Rmax constraint invalid\n");
    return false;
  }
  this->Rmax = Rmax;

  if (alpha <= 0.001 or alpha > 0.5){
    printf("MCChains::ERROR: alpha constraint invalid\n");
    return false;
  }
  this->alpha = alpha;
  
  return true;
}

bool MCChains::converged(){
  double R, CIt;
  double* pararray;
  double* totarray;
  int j,k,n,itot,cbase,upper,lower;
  size_t sortsize;
  int totlength=0;
  double CI,CIm=0;
  bool converged=false;

  for (i=0;i<nchains;i++)	
    totlength+=(chainlength[i]/2);
  totarray = new double[totlength];

  for (i=0;i<npars;i++){
    itot = 0;
    CIm = 0;
    for (j=0;j<nchains;j++){
      n = int(chainlength[j]/2);
      pararray = new double[n];
      cbase = j*chainwidth+i;
      for (k=n;k<chainlength[j];k++){
	totarray[itot] = chains[cbase][k];
	pararray[k-n] = chains[cbase][k];
	itot++;
      }
      //m chain math
      sortsize = size_t(n);
      qsort(pararray,sortsize,sizeof(double),dcomp);
      lower = (int)n*(alpha/2.0);
      upper = (int)n*(1.0-alpha/2.0);
      CI = abs((pararray[upper]-pararray[lower]));
      CIm += CI;
      printf("\nCI Chain %i: %f, L=%i, u=%i, b=%i\n",j,CI,n,upper,lower);
      delete[] pararray;
    }
    //m mean
    CIm = CIm/double(nchains);
    //t chain math
    sortsize = size_t(totlength);
    qsort(totarray,sortsize,sizeof(double),dcomp);
    lower = (int)totlength*(alpha/2.0);
    upper = (int)totlength*(1.0-alpha/2.0);
    CIt = abs((totarray[upper]-totarray[lower]));
    printf("CI mean: %f:\n",CIm);
    printf("CI Tot Chain: %f, L=%i, u=%i, b=%i\n",CIt,totlength,upper,lower);
    //r math
    R = CIt/CIm;
    printf("Param: %i, R: %f\n\n",i,R);
    converged = (converged and (R < Rmax));
  }
  delete[] totarray;
  
  return converged;
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
  char ctemp[2];
  string cnum;
  
  printf("Output Columns\n");
  for(int j=0;j<nchains;j++){
    sprintf(ctemp,"%i",j);
    cnum = string(ctemp);
    for (i=0;i<npars;i++){
      colnames[j*chainwidth+i] = parnames[i]+cnum;
    }
    colnames[j*chainwidth+npars] += cnum;
  }
  
  Table *newTable = pFits->addTable(hname,nruns,colnames,colform,colunits,AsciiTbl);
  
  for(i=0;i<allwidth;i++)
    printf("\t%s\n",colnames[i].c_str());
  newTable->column(colnames[i]).write(chains[i],1);
  return true;
}

MCChains::~MCChains(){
  delete[] bestpars;
  delete[] chainlength;
}
