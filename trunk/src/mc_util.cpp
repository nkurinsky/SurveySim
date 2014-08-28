/*
  Utility Classes Function Definitions for Monte Carlo Fitter
  Written by Noah Kurinsky, 10/31/13

  These classes are used for controlling acceptance probabilities of chain steps
  and tracking accepted moves for m chains
*/

#include "mc_util.h"

ParameterSettings::ParameterSettings(size_t nparams){
  min.resize(nparams);
  max.resize(nparams);
  sigma.resize(nparams);
  best.resize(nparams);
}

void ParameterSettings::set(short pnum, double Minimum, double Maximum, double standardDeviation, double bestValue){
  if(pnum < min.size()){
    min[pnum] = Minimum;
    max[pnum] = Maximum;
    sigma[pnum] = standardDeviation;
    best[pnum] = bestValue;
  }
  else
    printf("ParameterSettings::set called with incorrect pnum value %i\n",pnum);
}

int dcomp(const void * a,const void * b){
  
  if((*(double*)a) > (*(double*)b)) return 1;
  if((*(double*)a) < (*(double*)b)) return -1;
  else return 0;
}

MetropSampler::MetropSampler(int nchains, double maxTemp, double tempScale, double idealpct, double acpt_buf){
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
  tscale = tempScale;
  ideal_acceptance = idealpct;
  accept_buffer = acpt_buf;
  accept_ratio = 1.0/log(idealpct);
}

bool MetropSampler::accept(int chainnum, double trial){
  static double itemp,tester;
  tester=(trial-previous[chainnum])*tscale/temp;
  accepted = true;

  if(tester > 0){
    itemp=rng.flat(0.0,1.0);   //want uniform random number from 0-1
    accepted = (itemp < exp(-tester)) ? true : false;
  }
  
  iteration_total[chainnum]++;
  if (accepted){
    accept_total[chainnum]++;
    previous[chainnum] = trial;
  }
  
  recent.push_back(accepted);
  if(recent.size() > RECENT_NUM)
    recent.pop_front();
  
  trials.push_back(trial);
  if(trials.size() > RECENT_NUM)
    trials.pop_front();
  
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

double MetropSampler::acceptance_rate(){
  recent_num=0;
  for(unsigned long i=0;i<recent.size();i++)
    if(recent[i]) recent_num++; 
  
  return (double(recent_num)/double(recent.size()));
}

bool MetropSampler::anneal(){
  double rate = acceptance_rate();
  if(rate > (ideal_acceptance + accept_buffer)){
    temp--;
    return true;}
  else if (rate < (ideal_acceptance - accept_buffer)){
    temp++;
    return true;
  }
  else
    return false;
}

void MetropSampler::reset(){
  for (int i=0;i<nchains;i++){
    accept_total[i] = 0;
    iteration_total[i] = 0;
    previous[i] = 1.0E+4;
  }
  recent.clear();
}

MetropSampler::~MetropSampler(){
  delete[] accept_total;
  delete[] iteration_total;
  delete[] previous;
}


//chains implementation
MCChains::MCChains(int nchains, int npars, int nruns, int convstep){
  chainwidth = npars+1;
  allwidth = chainwidth*nchains;
  this->npars = npars;
  this->nruns = nruns;
  this->nchains = nchains;
  chains.resize(allwidth);
  for(i=0;i<allwidth;i++)
    chains[i].resize(nruns);
  rvals.resize(npars);
  convruns = nruns/convstep;
  for(i=0;i<npars;i++)
    rvals[i].resize(convruns);
  accepted.resize(nchains);
  for(i=0;i<nchains;i++)
    accepted[i].resize(nruns);
  bestpars = new double[npars];
  chainlength = new int[nchains];
  for (i=0;i<nchains;i++)
    chainlength[i]=0;
  chi_min = 1e10;
  //defaults to 95% CI
  Rmax = 1.05;
  alpha = 0.05;
}

bool MCChains::add_link(int chain, double pars[], double chisqr, bool accpt){
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
  accepted[chain][chainlength[chain]] = (accpt) ? 1.0 : 0.0;
  
  chainlength[chain]++;
  return true;
}

void MCChains::get_best_link(double pars[], double chisqr){
  for(int i=0;i<npars;i++)
    pars[i] = bestpars[i];
  chisqr = chi_min;
}

void MCChains::get_fit_results(double pars[], double sigma[]){
  vector<double> values;
  double mean, variance;
  unsigned long size;
  
  for(int i=0;i<npars;i++){
    variance = mean = 0.0;
    values.clear();
    for(int j=0;j<nchains;j++){
      size = chainlength[j];
      printf("\t%i %i %lu\n",i,j,size);
      for(unsigned long k = size/2; k < size; k++){
	values.push_back(chains[j*chainwidth+i][k]);
	mean += values.back();
      }
    }
    mean /= static_cast<double>(values.size());
    for(vector<double>::const_iterator val = values.begin(); val != values.end(); val++){
      variance += pow(*val - mean,2);
    }
    variance /= (static_cast<double>(values.size()) - 1);

    pars[i] = mean;
    sigma[i] = sqrt(variance);
  }
}

void MCChains::get_stdev(double sigma[]){
  double pars[npars];
  get_fit_results(pars,sigma);
}

bool MCChains::set_constraints(double Rmax, double alpha){
  bool retval=true;
  
  if (Rmax <= 1.0 or Rmax > 2.0){
    printf("MCChains::ERROR: Rmax constraint invalid\n");
    retval=false;
  }
  else
    this->Rmax = Rmax;
  
  if (alpha <= 0.001 or alpha > 0.5){
    printf("MCChains::ERROR: alpha constraint invalid\n");
    retval=false;
  }
  else
    this->alpha = alpha;
  
  return retval;
}

bool MCChains::converged(){
  static int call=0;
  double R, CIt;
  double* pararray;
  double* totarray;
  int j,k,n,itot,cbase,upper,lower;
  size_t sortsize;
  int totlength=0;
  double CI,CIm=0;
  bool converged=true;

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
      //printf("CI Chain %i: %f, L=%i, u=%i, b=%i\n",j,CI,n,upper,lower);
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
    //printf("CI mean: %f:\n",CIm);
    //printf("CI Tot Chain: %f, L=%i, u=%i, b=%i\n",CIt,totlength,upper,lower);
    //r math
    R = CIt/CIm;
    rvals[i][call] = R;
    printf("(%i) -  Param: %i, R: %f (%f)\n",call,i,R,Rmax);
    converged = (converged and (R < Rmax));
  }
  delete[] totarray;
  call++;

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
  
  std::slice domain = std::slice(0,chainlength[0],1);
  std::valarray<double> valarray_temp(chainlength[0]);
  
  int tablewidth = allwidth+nchains;
  std::vector<string> colnames(tablewidth,"CHISQ");
  std::vector<string> colunits(tablewidth,"-");
  std::vector<string> colform(tablewidth,"e13.5");
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
  for (int i=allwidth;i<tablewidth;i++){
    sprintf(ctemp,"%i",i-allwidth);
    cnum = string(ctemp);
    colnames[i] = "ACPT"+cnum;
  }

  printf("Making Table\n");
  Table *newTable = pFits->addTable(hname,chainlength[0],colnames,colform,colunits,AsciiTbl);
  printf("Made Table\n");

  int i; //I had to add this here, else wouldn't compile for some reason
  for(i=0;i<tablewidth;i++){
    if (i != 0)
      printf(", ");
    printf("%s",colnames[i].c_str());
    if(i < allwidth){
      valarray_temp = chains[i][domain];
      newTable->column(colnames[i]).write(valarray_temp,1);
    }
    else{
      valarray_temp = accepted[(i-allwidth)][domain];
      newTable->column(colnames[i]).write(valarray_temp,1);
    }
  }

  std::vector<string> colnames2(npars,"R");
  std::vector<string> colunits2(npars,"-");
  std::vector<string> colform2(npars,"f5.2");
  hname = "Convergence";
  
  printf("\nOutput Convergence Columns\n");
  for(i=0;i<npars;i++){
    sprintf(ctemp,"%i",i);
    cnum = string(ctemp);
    colnames2[i] += cnum;
  }
  
  newTable = pFits->addTable(hname,convruns,colnames2,colform2,colunits2,AsciiTbl);
  
  for(i=0;i<npars;i++){
    printf("%s,\t",colnames2[i].c_str());
    newTable->column(colnames2[i]).write(rvals[i],1);
  }
  
  printf("\n");
  return true;
}

MCChains::~MCChains(){
  delete[] bestpars;
  delete[] chainlength;
}

ResultChain::ResultChain(int num_arrays, int nresults){
  results.resize(num_arrays);
  chisqrs.reserve(nresults);
  for (int i=0;i<num_arrays;i++){
    results[i].reserve(nresults);
  }
}

bool ResultChain::add_link(valarray<double> arrays[], double chisqr){
  if(arrays == NULL){
    printf("NULL arrays\n");
    return false;}
  
  chisqrs.push_back(chisqr);
  for(i=0;i<results.size();i++){
    results[i].push_back(arrays[i]);
  }
  
  return true;
}

bool ResultChain::save(string filename, string resnames[]){
  using namespace CCfits;
  std::auto_ptr<FITS> pFits(0);

  try{
    pFits.reset(new FITS(filename,Write));
  }
  catch (CCfits::FITS::CantOpen){
    std::cerr << "Unknown Error Occurred, Can't save results" << endl;
    return false;
  }
  
  int tablewidth = int(results.size())+1;
  std::vector<string> colnames(tablewidth,"CHISQ");
  std::vector<string> colform(tablewidth,"1D");
  std::vector<string> colunit(tablewidth,"");
  
  colform[1] = "16D";
  colform[2] = "13D";
  colform[3] = "10D";

  for(unsigned long i=0;i<results.size();i++)
    colnames[i+1]=resnames[i];

  string hname("Results Extension");
  int nrows(chisqrs.size());

  printf("Making Table: %i\n",nrows);
  Table *newTable = pFits->addTable(hname,nrows,colnames,colform,colunit);
  printf("Made Table\n");

  newTable->column(colnames[0]).write(chisqrs,1);
  printf("%s ",colnames[0].c_str());
  for(i=0;i<results.size();i++){
    printf("%s ",colnames[i+1].c_str());
    newTable->column(colnames[i+1]).writeArrays(results[i],1);
  }

  printf("\n");
  return true;
}