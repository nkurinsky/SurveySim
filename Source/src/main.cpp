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

using namespace std; 

int main(int argc,char** argv){
  
  printf("\nMCMC Fitting Code Compiled: %s\n\n",__DATE__);
  
  Configuration q(argc,argv);
  q.print();
  RandomNumberGenerator rng;
  
  //general variable declarations
  bool accept,saved;
  unsigned long i,m;
  unsigned long pi;
  vector<int>::const_iterator p;
  double temp, trial;

  //this array holds phi0,L0,alpha,beta,p, and q
  double lpars[LUMPARS];
  double lmin[LUMPARS];
  double lmax[LUMPARS];
  double lsigma[LUMPARS];
  q.LFparameters(lpars);
  q.LFparameters(lmin,q.min);
  q.LFparameters(lmax,q.max);
  q.LFparameters(lsigma,q.step);
  double CE = q.colorEvolution[q.value];
  double CEmin = q.colorEvolution[q.min];
  double CEmax = q.colorEvolution[q.max];
  double CEsigma = q.colorEvolution[q.step];
  
  //counts
  string countnames[]= {"dnds250","dnds350","dnds500"};

  /*
    note that the dp values here are the widths of the "proposal distribution"
    the smaller they are the slower will converge onto the right answer
    the bigger they are, the less well sampled the probability distribution 
    will be and hence the less accurate the final answer
    needs to find the "goldilocks" zone here but that's very problem specific 
    so requires trial and error as we actually start fitting data. 
  */

  lumfunct lf;
  lf.set_params(lpars);

  //initialize simulator
  simulator survey(q.filterfile,q.obsfile,q.sedfile); 
  survey.set_color_exp(CE); //color evolution
  survey.set_lumfunct(&lf);
  
  //NON GENERAL COUNTS
  //the flux array is logarithmic in steps of 0.3dex 
  q.ns=8;
  survey.set_size(q.areaSteradian(),q.dz,q.zmin,q.nz,q.ns);
  int bnum[]={16,13,10};
  products output(q.nz,bnum);

  /*
    Loop trough the mcmc runs
    At each step, from some initial paramer value "p" call on the simulator.cpp 
    program to evaluate the chi2 for that particular color-color plot. 
    Then use the metrop algorithm (below) to decide whether or not to keep a particular guess
    The chain results (including all accepted guesses) will provide us with the posterior 
    probability distributions for the fitted parameters.
  */
  
  MCChains mcchain(q.nchain,q.nparams,q.runs,q.conv_step);
  ResultChain counts(3,q.nchain*q.runs);
  mcchain.set_constraints(q.rmax,q.a_ci);
  MetropSampler metrop(q.nchain,q.tmax,q.idealpct,q.annrng);

  //mcmc variables and arrays
  double chi_min=1.0E+4; 
  double ptemp[q.nchain][q.nparams];
  double pbest[q.nparams];
  double pcurrent[q.nchain][q.nparams];
  

  //initialize first chain to initial parameters
  for(pi=0, p = q.param_inds.begin(); p != q.param_inds.end(); ++p,++pi)
    ptemp[0][pi] = pcurrent[0][pi] = lpars[*p];
  if(q.vary_cexp)
    ptemp[0][q.cind] = pcurrent[0][q.cind] = CE;
  
  for(m=1;m<q.nchain;m++){ 
    for(pi=0, p = q.param_inds.begin(); p != q.param_inds.end(); ++p,++pi)
      ptemp[m][pi] = pcurrent[m][pi] = rng.flat(lmin[*p],lmax[*p]);
    if(q.vary_cexp)
      ptemp[m][q.cind] = pcurrent[m][q.cind] = rng.flat(CEmin,CEmax);
  }
  
  if (q.oprint)
    printf("\n ---------- Beginning MC Burn-in Phase ---------- \n");
  else
    printf("\nBeginning MC Burn-in Phase...\n");

  for (i=0;i<q.burn_num;i++){
    for (m=0;m<q.nchain;m++){
      for(pi=0, p = q.param_inds.begin(); p != q.param_inds.end(); ++p,++pi){
	temp=rng.gaussian(pcurrent[m][pi],lsigma[*p]);
	if( ( temp >= lmin[*p] ) and ( temp <= lmin[*p] ) ) {
	  ptemp[m][pi]=temp;
	  lpars[*p] = temp;
	}
      }
      if(q.vary_cexp){
	temp=rng.gaussian(pcurrent[m][q.cind],CEsigma);
	if( ( temp >= CEmin ) and ( temp <= CEmax ) ) {
	  ptemp[m][q.cind]=temp;
	  survey.set_color_exp(temp);
	}
      }

      if (q.oprint) 
	printf("%lu %lu - %lf %lf %lf %lf %lf %lf %lf %lf : ",(i+1),(m+1),lpars[0],lpars[1],lpars[2],lpars[3],lpars[4],lpars[5],lpars[6],(q.vary_cexp ? ptemp[m][q.cind] : CE));
      lf.set_params(lpars);
      output=survey.simulate();
      trial=output.chisqr;
      if (q.oprint)
	printf("Iteration Chi-Square: %lf",trial);
      
      if(trial < chi_min){
	if (q.oprint) 
	  printf(" -- Current Best Trial");
	chi_min=trial;
	for(pi=0;pi<q.nparams;pi++)
	  pbest[pi]=ptemp[m][pi];
      }
      if (q.oprint)
	printf("\n");

      if(metrop.accept(m,trial)){
	for(pi=0;pi<q.nparams;pi++)
	  pcurrent[m][pi]=ptemp[m][pi];
      }      
    }
    if(((i+1) % q.burn_step) == 0)
      if(not metrop.anneal())
	i = q.runs;
  }
  
  for(m=0;m<q.nchain;m++)
    for(pi=0;pi<q.nparams;pi++)
      ptemp[m][pi] = pcurrent[m][pi] = pbest[pi];
  
  metrop.reset();
  
  if (q.oprint)
    printf("\n\n ---------------- Fitting Start ---------------- \n Total Run Number: %ld\n\n",q.runs);
  else
    printf("\nBeginning Fitting, Run Maximum: %ld\n",q.runs);
  
  for (i=0;i<q.runs;i++){
    for (m=0;m<q.nchain;m++){
      for(pi=0, p = q.param_inds.begin(); p != q.param_inds.end(); ++p,++pi){
	temp = rng.gaussian(pcurrent[m][pi],lsigma[*p]);
	if( ( temp >= lmin[*p] ) and ( temp <= lmax[*p] ) ) {
	  ptemp[m][pi]=temp;
	  lpars[*p] = temp;
	}
      }

      if(q.vary_cexp){
	temp = rng.gaussian(pcurrent[m][q.cind],CEsigma);
	if( (temp >= CEmin ) and ( temp <= CEmax ) ) {
	  ptemp[m][q.cind] = temp;
	  survey.set_color_exp(temp);
	}
      }
      
      if (q.oprint){
	printf("%lu %lu -",(i+1),(m+1));
	for(pi=0;pi<q.nparams;pi++)
	  printf(" %lf",ptemp[m][pi]);
	printf(" : ");
      }
      lf.set_params(lpars);
      output=survey.simulate();
      trial=output.chisqr;
      if (q.oprint)
	printf("Model Chi-Square: %lf",trial);
      
      accept = metrop.accept(m,trial);
      mcchain.add_link(m,ptemp[m],trial,accept);
      counts.add_link(output.dnds,trial);

      if(accept){
	if (q.oprint) 
	  printf(" -- Accepted\n");
	for(pi=0;pi<q.nparams;pi++)
	  pcurrent[m][pi]=ptemp[m][pi];
      }
      else if (q.oprint)
	printf(" -- Rejected\n");
    }
    
    if(((i+1) % q.conv_step) == 0){
      printf("Checking Convergence\n");
      if (i < q.burn_num)
	metrop.anneal();
      if(mcchain.converged()){
	printf("Chains Converged!\n");
	i = q.runs;}
      else
	printf("Chains Have Not Converged\n");
    } 
  }
  
  mcchain.get_best_link(pbest,chi_min);

  for(pi=0;pi<q.param_inds.size();pi++)
    lpars[q.param_inds[pi]] = pbest[pi];
  lf.set_params(lpars);
  if(q.vary_cexp)
    survey.set_color_exp(pbest[q.cind]);

  printf("\nRe-Running Best Fit...\n");
  output=survey.simulate();
  printf("Model chi2: %lf\n",output.chisqr);
  printf("Acceptance Rate: %lf%%\n",metrop.mean_acceptance());
  
  string pnames[] = {"PHI0","L0","ALPHA","BETA","P","Q","ZCUT"}; 
  unique_ptr<string[]> parnames (new string[q.nparams]);
  for(pi=0;pi<q.param_inds.size();pi++)
    parnames[pi] = pnames[q.param_inds[pi]];
  if(q.vary_cexp)
    parnames[q.cind] = "CEXP";
  
  saved = survey.save(q.outfile);
  saved &= mcchain.save(q.outfile,parnames.get());
  saved &= counts.save(q.outfile,countnames);
  if(saved)
    printf("Save Successful\n");
  else
    printf("Save Failed\n");
  
  printf("\nFitting Complete\n\n");
  
  if(saved)
    return 0;
  else
    return 1;
}
