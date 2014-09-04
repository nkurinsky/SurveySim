//****************************************************************************
// written by Anna Sajina and Noah Kurinsky, 2012
// updated by Noah Kurinsky 10/13
// The purpose of this program is to read-in the parameters set in the idl 
// wrapper (simulation.pro) and to use MCMC to determine the best-fit set of 
// parameters. 
// Specifically it runs through a number of iterations (runs) in each stage 
// calling on simulator.cpp to compute the chi2 for a given set of parameters
// then the metrop algorithm is invoked to decide whether or not a given trial 
// run is accepted. If so it is added to a chain of mcmc values.
//****************************************************************************

#include "simulator.h"
#include "functions.h"
#include "mc_util.h"
#include <stdio.h>

//Array size definitions 
#define BANDS 3

#define VERBOSE( ARG ) \
  if(q.oprint){	       \
    ARG;	       \
  }

#define TERSE( ARG )   \
  if(!q.oprint){       \
    ARG;	       \
  }

using namespace std; 

int main(int argc,char** argv){
  
  printf("\nMCMC Fitting Code Compiled: %s\n\n",__DATE__);
  
  Configuration q(argc,argv);
  VERBOSE(q.print());
  RandomNumberGenerator rng;
  
  //general variable declarations
  bool accept,saved;
  unsigned long i,m;
  unsigned long pi;
  double trial;

  //this array holds phi0,L0,alpha,beta,p, and q
  double lpars[LUMPARS];
  double lmin[LUMPARS];
  double lmax[LUMPARS];
  q.LFparameters(lpars);
  q.LFparameters(lmin,q.min);
  q.LFparameters(lmax,q.max);
  double CE = q.colorEvolution[q.value];
  double CEmin = q.colorEvolution[q.min];
  double CEmax = q.colorEvolution[q.max];
  
  //counts
  string countnames[]= {"dnds250","dnds350","dnds500"};
  
  lumfunct lf;
  lf.set_params(lpars);

  //initialize simulator
  simulator survey(q.filterfile,q.obsfile,q.sedfile); 
  survey.set_color_exp(CE); //color evolution
  survey.set_lumfunct(&lf);
  
  //NON GENERAL COUNTS
  //the flux array is logarithmic in steps of 0.3dex 
  survey.set_size(q.areaSteradian(),q.dz,q.zmin,q.nz);
  products output;
  
  MCChains mcchain(q.nchain,q.nparams,q.runs,q.conv_step);
  MCChains burnchain(q.nchain,q.nparams,q.runs,q.conv_step);
  mcchain.set_constraints(q.rmax,q.a_ci);
  burnchain.set_constraints(q.rmax,q.a_ci);

  ResultChain counts(3,q.nchain*q.runs);
 
  MetropSampler metrop(q.nchain,q.tmax,q.tscale,q.idealpct,q.annrng);

  //mcmc variables and arrays
  double chi_min=1.0E+4; 
  double ptemp[q.nchain][q.nparams];
  double pcurrent[q.nchain][q.nparams];
  double *setvars[q.nparams];
  ParameterSettings pset(q.nparams);

  VERBOSE(printf("Initializing Chains\n"));
  vector<int>::const_iterator p;
  for(pi=0, p = q.param_inds.begin(); p != q.param_inds.end(); ++p,++pi){
    setvars[pi] = &(lpars[*p]);
    pcurrent[0][pi] = *(setvars[pi]);
    pset.set(pi,lmin[*p],lmax[*p],(lmax[*p] - lmin[*p])/6.0,lpars[*p]);
  }
  if(q.vary_cexp){
    setvars[q.cind] = &CE;
    pcurrent[0][q.cind] = *(setvars[pi]);
    pset.set(q.cind,CEmin,CEmax,(CEmax - CEmin)/6.0,CE);
  }

  for(m=1;m<q.nchain;m++){ 
    for(pi=0;pi<q.nparams;pi++)
      pcurrent[m][pi] = rng.flat(pset.min[pi],pset.max[pi]);
  }
  
  VERBOSE(printf("\n ---------- Beginning MC Burn-in Phase ---------- \n"));
  TERSE(printf("\nBeginning MC Burn-in Phase...\n"));

  for (i=0;i<q.burn_num;i++){
    for (m=0;m<q.nchain;m++){
      
      for(pi = 0;pi<q.nparams;pi++)
	*(setvars[pi]) = ptemp[m][pi] = rng.gaussian(pcurrent[m][pi],pset.sigma[pi],pset.min[pi],pset.max[pi]);
      
      VERBOSE(printf("%lu %lu -",(i+1),(m+1)));
      for(pi=0;pi<LUMPARS;pi++)
	VERBOSE(printf(" %5.2lf",lpars[pi]));
      VERBOSE(printf(" %5.2lf : ",(q.vary_cexp ? ptemp[m][q.cind] : CE)));
      
      lf.set_params(lpars);
      survey.set_color_exp(CE);
      
      output=survey.simulate();
      trial=output.chisqr;
      
      VERBOSE(printf("Iteration Chi-Square: %lf",trial));      
      if(trial < chi_min){
	VERBOSE(printf(" -- Current Best Trial"));
	chi_min=trial;
	for(pi=0;pi<q.nparams;pi++)
	  pset.best[pi]=ptemp[m][pi];
      }
      VERBOSE(printf("\n"));

      accept = metrop.accept(m,trial);
      if(accept)
	for(pi=0;pi<q.nparams;pi++)
	  pcurrent[m][pi]=ptemp[m][pi];

      burnchain.add_link(m,ptemp[m],trial,accept);
    }
    
    if(((i+1) % q.burn_step) == 0){
      printf("\nAcceptance: %5.1lf%%\n",metrop.acceptance_rate()*100.0);
      if(not metrop.anneal())
	i = q.runs;
      burnchain.get_stdev(pset.sigma.data());
    }
  }
  
  for(m=0;m<q.nchain;m++)
    for(pi=0;pi<q.nparams;pi++)
      pcurrent[m][pi] = pset.best[pi];
  
  metrop.reset();
  
  VERBOSE(printf("\n\n ---------------- Fitting Start ---------------- \n Total Run Number: %ld\n\n",q.runs));
  TERSE(printf("\nBeginning Fitting, Run Maximum: %ld\n",q.runs));
  
  for (i=0;i<q.runs;i++){
    for (m=0;m<q.nchain;m++){
      
      for(pi = 0;pi<q.nparams;pi++)
	*(setvars[pi]) = ptemp[m][pi] = rng.gaussian(pcurrent[m][pi],pset.sigma[pi],pset.min[pi],pset.max[pi]);
      
      VERBOSE(printf("%lu %lu -",(i+1),(m+1)));
      for(pi=0;pi<q.nparams;pi++)
	VERBOSE(printf(" %lf",ptemp[m][pi]));
      VERBOSE(printf(" : "));
      
      lf.set_params(lpars);
      survey.set_color_exp(CE);

      output=survey.simulate();
      trial=output.chisqr;
      VERBOSE(printf("Model Chi-Square: %lf",trial));
      
      accept = metrop.accept(m,trial);

      mcchain.add_link(m,ptemp[m],trial,accept);
      counts.add_link(output.dnds,trial);

      if(accept){
	VERBOSE(printf(" -- Accepted\n"));
	for(pi=0;pi<q.nparams;pi++)
	  pcurrent[m][pi]=ptemp[m][pi];
      }
      else
	VERBOSE(printf(" -- Rejected\n"));
    }
    
    if(((i+1) % q.conv_step) == 0){
      printf("Checking Convergence\n");
      if (i < q.burn_num)
	metrop.anneal();
      mcchain.get_stdev(pset.sigma.data());
      if(mcchain.converged()){
	printf("Chains Converged!\n");
	i = q.runs;}
      else
	printf("Chains Have Not Converged\n");
    } 
  }
  
  mcchain.get_best_link(pset.best.data(),chi_min);

  string pnames[] = {"PHI0","L0","ALPHA","BETA","P","Q","ZCUT"};
  printf("Best fit:\n");
  for(pi=0;pi<q.param_inds.size();pi++){
    lpars[q.param_inds[pi]] = pset.best[pi];
    printf("%s : %lf +\\- %lf\n",pnames[q.param_inds[pi]].c_str(),pset.best[pi],pset.sigma[pi]);
  }
  lf.set_params(lpars);
  if(q.vary_cexp){
    survey.set_color_exp(pset.best[q.cind]);
    printf("CEXP : %lf +\\- %lf\n",pset.best[pi],pset.sigma[pi]);
  }  

  printf("\nRe-Running Best Fit...\n");
  output=survey.simulate();
  printf("Model chi2: %lf (%lf)\n",output.chisqr,chi_min);
  printf("Acceptance Rate: %lf%%\n",metrop.mean_acceptance()*100);

  unique_ptr<string[]> parnames (new string[q.nparams]);
  
  for(pi=0, p= q.param_inds.begin(); p != q.param_inds.end();pi++,p++)
    parnames[pi] = pnames[*p];
  if(q.vary_cexp)
    parnames[q.cind] = "CEXP";
  
  VERBOSE(printf("Saving Survey\n"));
  saved = survey.save(q.outfile);
  VERBOSE(printf("Saving Chains\n"));
  saved &= mcchain.save(q.outfile,parnames.get());
  VERBOSE(printf("Saving Counts\n"));
  saved &= counts.save(q.outfile,countnames);

  saved ? printf("Save Successful") : printf("Save Failed");
  printf("\nFitting Complete\n\n");
  
  return saved ? 0 : 1;
}


