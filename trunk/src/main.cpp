//****************************************************************************
// SurveySim
//.............................................................................
// written by Noah Kurinsky and Anna Sajina
// last modified Anna Sajina November 1st, 2014
// The purpose of this program is to use MCMC to determine the luminosity function evolution and redshift distribution based on infrared survey data -- the details of the input observations as well as default model parameters are passed to this program via a series of fits files (model.fits, observations.fits, sed_templates.fits) which are set using either idl (SurveySim.pro) or python (SurveySim.py). 
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
  printf("\n");

  Configuration q(argc,argv);

  VERBOSE(printf("Fitting Code Compiled: %s\n\n",__DATE__));
  VERBOSE(q.print());
  RandomNumberGenerator rng;

  printf("Axis 1: %s\nAxis 2: %s\n",get_axis_type(q.axes[0]).c_str(),get_axis_type(q.axes[1]).c_str());

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
  
  string pnames[] = {"PHI0","L0","ALPHA","BETA","P","Q","ZCUT"};
  string countnames[]= {"dnds1","dnds2","dnds3"};
  
  lumfunct lf;
  lf.set_params(lpars);

  //initialize simulator
  simulator survey(q); 
  survey.set_color_exp(CE); //color evolution
  survey.set_lumfunct(&lf);
  products output;

  MCChains mcchain(q.nchain,q.nparams,q.runs,q.conv_step);
  MCChains burnchain(q.nchain,q.nparams,q.runs,q.conv_step);
  mcchain.set_constraints(q.rmax,q.a_ci);
  burnchain.set_constraints(q.rmax,q.a_ci);

  ResultChain counts(3,q.nchain*q.runs);
  ResultChain final_counts(3,q.nsim);
 
  MetropSampler metrop(q.nchain,q.tmax,q.tscale,q.idealpct,q.annrng);

  //mcmc variables and arrays
  double chi_min=1.0E+4; 
  vector<int>::const_iterator p;

  double *setvars[q.nparams];
  ParameterSettings pset(q.nparams);
  
  if(q.nparams > 0){
    double ptemp[q.nchain][q.nparams];
    vector<double > means(q.nparams);
    double pcurrent[q.nchain][q.nparams];
    
    VERBOSE(printf("Initializing Chains\n"));
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

	vector<double> results;
	for(pi=0;pi<q.nparams;pi++)
	  means[pi] = pcurrent[m][pi];
	rng.gaussian_mv(means,pset.covar,pset.min,pset.max,results);
	for(pi = 0;pi<q.nparams;pi++){
	  //*(setvars[pi]) = ptemp[m][pi] = rng.gaussian(pcurrent[m][pi],pset.sigma[pi],pset.min[pi],pset.max[pi]);
	  *(setvars[pi]) = ptemp[m][pi] = results[pi];
	}
	
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
	printf("\nAcceptance: %5.1lf%% (T=%8.2e)\n",metrop.acceptance_rate()*100.0,metrop.temperature());
	printf("Covariance Matrix:\n");
	for(int pi=0;pi<pset.covar.size();pi++){
	  printf("\t[");
	  for(int pj=0;pj<pset.covar[pi].size();pj++)
	    printf("%6.2lf ",pset.covar[pi][pj]);
	  printf("]\n");
	}
	if(not metrop.anneal())
	  i = q.runs;
	burnchain.get_stdev(pset.sigma.data());
	burnchain.get_covariance(pset.covar);
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

	vector<double> results;
	for(pi = 0;pi<q.nparams;pi++)
	  means[pi] = pcurrent[m][pi];
	rng.gaussian_mv(means,pset.covar,pset.min,pset.max,results);
	for(pi = 0;pi<q.nparams;pi++){
	  //*(setvars[pi]) = ptemp[m][pi] = rng.gaussian(pcurrent[m][pi],pset.sigma[pi],pset.min[pi],pset.max[pi]);
	  *(setvars[pi]) = ptemp[m][pi] = results[pi];
	}
		
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
	mcchain.get_covariance(pset.covar);
	if(mcchain.converged()){
	  printf("Chains Converged!\n");
	  i = q.runs;}
	else
	  printf("Chains Have Not Converged\n");
      } 
    }
    
    mcchain.get_best_link(pset.best.data(),chi_min);
    
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
    printf("\nCovariance Matrix:\n");
    for(int pi=0;pi<pset.covar.size();pi++){
      printf("\t[");
      for(int pj=0;pj<pset.covar[pi].size();pj++)
	printf("%6.2lf ",pset.covar[pi][pj]);
      printf("]\n");
    }
    
    printf("Final Acceptance Rate: %lf%%\n",metrop.mean_acceptance()*100);
  }

  printf("\nRunning Best Fit...\n");
  output=survey.simulate();
  double tchi_min = output.chisqr;
  VERBOSE(printf("Saving Initial Survey\n"));
  saved = survey.save(q.outfile);

  vector<double> results,means(q.nparams);
  for(pi = 0;pi<q.nparams;pi++)
    means[pi] = pset.best[pi];
  
  for(unsigned long i=1;i<q.nsim;i++){

    //vary according to final results

    rng.gaussian_mv(means,pset.covar,pset.min,pset.max,results);
    for(pi = 0;pi<q.nparams;pi++){
      //*(setvars[pi]) = rng.gaussian(pset.best[pi],pset.sigma[pi],pset.min[pi],pset.max[pi]);
      *(setvars[pi]) = results[pi];
    }
    
    lf.set_params(lpars);
    survey.set_color_exp(CE);
    
    output=survey.simulate();
    if(output.chisqr < tchi_min){
      tchi_min = output.chisqr;
      printf("Saving new minimum (%lu/%lu)...",i,q.nsim);
      saved = survey.save(q.outfile);
    }
    final_counts.add_link(output.dnds,output.chisqr);
  }
  printf("Saved Chi2: %lf (%lf)\n",tchi_min,chi_min);
  
  unique_ptr<string[]> parnames (new string[q.nparams]);
  
  for(pi=0, p= q.param_inds.begin(); p != q.param_inds.end();pi++,p++)
    parnames[pi] = pnames[*p];
  if(q.vary_cexp)
    parnames[q.cind] = "CEXP";
  
  VERBOSE(printf("Saving Chains\n"));
  saved &= mcchain.save(q.outfile,parnames.get());
  //VERBOSE(printf("Saving Counts\n"));
  saved &= counts.save(q.outfile,countnames);
  VERBOSE(printf("Saving Final Counts\n"));
  saved &= final_counts.save(q.outfile,countnames,"Counts Monte Carlo");

  saved ? printf("Save Successful") : printf("Save Failed");
  printf("\nFitting Complete\n\n");
  
  return saved ? 0 : 1;
}


