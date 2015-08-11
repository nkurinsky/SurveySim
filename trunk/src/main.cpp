//****************************************************************************
// SurveySim
//.............................................................................
// written by Noah Kurinsky and Anna Sajina
// The purpose of this program is to use MCMC to determine the luminosity function evolution and redshift distribution based on infrared survey data -- the details of the input observations as well as default model parameters are passed to this program via a series of fits files (model.fits, observations.fits, templates.fits) which are set using either idl (SurveySim.pro) or python (SurveySim.py). 
//****************************************************************************

#include "simulator.h"
#include "functions.h"
#include "mc_util.h"
#include <stdio.h>

//Array size definitions 
#define BANDS 3

using namespace std; 

int main(int argc,char** argv){
  Configuration q(argc,argv);

  int logflag=q.oprint;

  LOG_DEBUG(printf("Fitting Code Compiled: %s\n\n",__DATE__));
  RandomNumberGenerator rng;

  LOG_DEBUG(printf("Axis 1: %s\nAxis 2: %s\n",get_axis_type(q.axes[0]).c_str(),get_axis_type(q.axes[1]).c_str()));

  if(logflag >= 2) q.print();

  //general variable declarations
  bool accept;
  bool saved=false;
  unsigned long i,m;
  unsigned long pi;
  double trial;

  //this array holds phi0,L0,alpha,beta,p, and q
  double lpars[LUMPARS];
  double lmin[LUMPARS];
  double lmax[LUMPARS];
  double parfixed[LUMPARS];

  //printf(q.LFparameters)
  q.LFparameters(lpars);
  q.LFparameters(lmin,q.min);
  q.LFparameters(lmax,q.max);
 
  double CE = q.colorEvolution[q.value];
  double CEmin = q.colorEvolution[q.min];
  double CEmax = q.colorEvolution[q.max];
  double ZBC = q.colorZCut[q.value];
  double ZBCmin = q.colorZCut[q.min];
  double ZBCmax = q.colorZCut[q.max];
  q.LFparameters(parfixed,q.fixed);

  string pnames[] = {"PHI0","L0","ALPHA","BETA","P","Q","P2","Q2","ZBP","ZBQ"};  
  string countnames[]= {"dnds1","dnds2","dnds3"};
  
  lumfunct lf(q.lfDist);
  lf.set_params(lpars);

  //initialize simulator
  simulator survey(q); 
  survey.set_color_exp(CE,ZBC); //color evolution
  survey.set_fagn_pars(lpars);
  survey.set_lumfunct(&lf);
  products output;
  ResultChain final_counts(3,q.nsim);

  if(q.simflag){ //observations not provided
    LOG_CRITICAL(printf("Simulating a survey only \n"));
    //run number of simulations to get counts range

    LOG_CRITICAL(printf("Beginning Simulation Loop (%lu runs)\n",q.nsim));
    for(unsigned long simi=1;simi<q.nsim;simi++){      
      output=survey.simulate();
      final_counts.add_link(output.dnds,output.chisqr);
      if((simi % 20) == 0){
	LOG_CRITICAL(printf("Current Iteration: %lu\n",simi));
      }
    }
    
    LOG_DEBUG(printf("\nSaving Output\n"));
    saved = survey.save(q.outfile);
    LOG_DEBUG(printf("Saving Counts\n"));
    saved &= final_counts.save(q.outfile,countnames,"Simulated Counts");
  }
  else{    

    MCChains mcchain(q.nchain,q.nparams,q.runs,q.conv_step);
    MCChains burnchain(q.nchain,q.nparams,q.runs,q.conv_step);
    mcchain.set_constraints(q.rmax,q.a_ci);
    burnchain.set_constraints(q.rmax,q.a_ci);
    
    //mcmc variables and arrays
    double chi_min=1.0E+4; 
    double tchi_min=chi_min;
    vector<int>::const_iterator p;
    
    double *setvars[q.nparams];
    ParameterSettings pset(q.nparams);
    
    unique_ptr<string[]> parnames (new string[q.nparams]);
    vector<double> results,means(q.nparams);
    
    ResultChain counts(3,q.nchain*q.runs);
    MetropSampler metrop(q.nchain,q.temp,q.learningRate,q.idealpct,q.annrng,q.oprint);    

    if(q.nparams > 0){ //observations provided, not all parameters fixed
      double ptemp[q.nchain][q.nparams];
      vector<double > means(q.nparams);
      double pcurrent[q.nchain][q.nparams];
      
      LOG_DEBUG(printf("Initializing Chains\n"));
      for(pi=0, p = q.param_inds.begin(); p != q.param_inds.end(); ++p,++pi){
	setvars[pi] = &(lpars[*p]);
	pcurrent[0][pi] = *(setvars[pi]);
	pset.set(pi,lmin[*p],lmax[*p],(lmax[*p] - lmin[*p])/6.0,lpars[*p]);
      }
      if(q.vary_cexp){
	setvars[q.cind] = &CE;
	pcurrent[0][q.cind] = *(setvars[pi]);
	pset.set(q.cind,CEmin,CEmax,(CEmax - CEmin)/6.0,CE);
	pi++;
      }
      if(q.vary_zbc){
	setvars[q.zbcind] = &ZBC;
	pcurrent[0][q.zbcind] = *(setvars[pi]);
	pset.set(q.zbcind,ZBCmin,ZBCmax,(ZBCmax - ZBCmin)/6.0,ZBC);
      }
      
      for(m=1;m<q.nchain;m++){ 
	for(pi=0;pi<q.nparams;pi++)
	  pcurrent[m][pi] = rng.flat(pset.min[pi],pset.max[pi]);
      }
      
      LOG_INFO(printf("\n ---------- Beginning MC Burn-in Phase ---------- \n"));
      i=0;
      for (i=0;i<q.burn_num;i++){
	for (m=0;m<q.nchain;m++){
	  
	  vector<double> results;
	  for(pi=0;pi<q.nparams;pi++)
	    means[pi] = pcurrent[m][pi];
	  rng.gaussian_mv(means,pset.covar,pset.min,pset.max,results);
	  for(pi = 0;pi<q.nparams;pi++){
	    *(setvars[pi]) = ptemp[m][pi] = results[pi];
	  }
	  
	  LOG_DEBUG(printf("%lu %lu -",(i+1),(m+1)));
	  for(pi=0;pi<LUMPARS;pi++)
	    LOG_DEBUG(printf(" %5.2lf",lpars[pi]));
	  LOG_DEBUG(printf(" %5.2lf : ",(q.vary_cexp ? ptemp[m][q.cind] : CE)));
	  LOG_DEBUG(printf(" %5.2lf : ",(q.vary_zbc ? ptemp[m][q.zbcind] : ZBC)));
	  
	  lf.set_params(lpars);
	  survey.set_color_exp(CE,ZBC);
	  survey.set_fagn_pars(lpars);
	  
	  output=survey.simulate();
	  
	  trial=output.chisqr;
	  
	  LOG_DEBUG(printf("Iteration Chi-Square: %lf",trial));      
	  if(trial < chi_min){
	    LOG_DEBUG(printf(" -- Current Best Trial"));
	    chi_min=trial;
	    for(pi=0;pi<q.nparams;pi++)
	      pset.best[pi]=ptemp[m][pi];
	  }
	  LOG_DEBUG(printf("\n"));
	  
	  accept = metrop.accept(m,trial);
	  if(accept)
	    for(pi=0;pi<q.nparams;pi++)
	      pcurrent[m][pi]=ptemp[m][pi];
	  
	  burnchain.add_link(m,ptemp[m],trial,accept);
	}
	
	if(((i+1) % q.burn_step) == 0){
	  LOG_INFO(printf("\nAcceptance: %5.1lf%% (T=%8.2e)\n",metrop.acceptance_rate()*100.0,metrop.temperature()));
	  LOG_INFO(printf("Covariance Matrix:\n"));
	  for(int pi=0;pi<pset.covar.size();pi++){
	    LOG_INFO(printf("\t["));
	    for(int pj=0;pj<pset.covar[pi].size();pj++)
	      LOG_INFO(printf("%6.2lf ",pset.covar[pi][pj]));
	    LOG_INFO(printf("]\n"));
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
      
      LOG_INFO(printf("\n\n ---------------- Fitting Start ---------------- \n Total Run Number: %ld\n\n",q.runs));
      
      for (i=0;i<q.runs;i++){
	for (m=0;m<q.nchain;m++){
	  
	  vector<double> results;
	  for(pi = 0;pi<q.nparams;pi++)
	    means[pi] = pcurrent[m][pi];
	  rng.gaussian_mv(means,pset.covar,pset.min,pset.max,results);
	  for(pi = 0;pi<q.nparams;pi++){
	    *(setvars[pi]) = ptemp[m][pi] = results[pi];
	  }
	  
	  LOG_DEBUG(printf("%lu %lu -",(i+1),(m+1)));
	  for(pi=0;pi<q.nparams;pi++)
	    LOG_DEBUG(printf(" %lf",ptemp[m][pi]));
	    LOG_DEBUG(printf(" : "));
	  
	  lf.set_params(lpars);
	  survey.set_color_exp(CE,ZBC);
	  survey.set_fagn_pars(lpars);
	  
	  output=survey.simulate();
	  trial=output.chisqr;
	  LOG_DEBUG(printf("Model Chi-Square: %lf",trial));
	  
	  accept = metrop.accept(m,trial);
	  
	  mcchain.add_link(m,ptemp[m],trial,accept);
	  counts.add_link(output.dnds,trial);
	  
	  if(accept){
	    LOG_DEBUG(printf(" -- Accepted\n"));
	    for(pi=0;pi<q.nparams;pi++)
	      pcurrent[m][pi]=ptemp[m][pi];
	  }
	  else
	    LOG_DEBUG(printf(" -- Rejected\n"));
	}
	
	if(((i+1) % q.conv_step) == 0){
	  LOG_INFO(printf("Checking Convergence\n"));
	  if (i < q.burn_num)
	    metrop.anneal();
	  mcchain.get_stdev(pset.sigma.data());
	  mcchain.get_covariance(pset.covar);
	  if(mcchain.converged()){
	    LOG_INFO(printf("Chains Converged!\n"));
	    i = q.runs;}
	  else
	    LOG_INFO(printf("Chains Have Not Converged\n"));
	} 
      }
      
      mcchain.get_best_link(pset.best.data(),chi_min);
      
      LOG_INFO(printf("Best fit:\n"));
      for(pi=0;pi<q.param_inds.size();pi++){
	lpars[q.param_inds[pi]] = pset.best[pi];
	LOG_INFO(printf("%s : %lf +\\- %lf\n",pnames[q.param_inds[pi]].c_str(),pset.best[pi],pset.sigma[pi]));
      }
      lf.set_params(lpars);
      if(q.vary_cexp){
	CE=pset.best[q.cind];
	LOG_DEBUG(printf("CEXP : %lf +\\- %lf\n",pset.best[pi],pset.sigma[pi]));
	pi++;
      }  
      if(q.vary_zbc){
	ZBC=pset.best[q.zbcind];
	LOG_DEBUG(printf("ZBC : %lf +\\- %lf\n",pset.best[pi],pset.sigma[pi]));
      }
      survey.set_color_exp(CE,ZBC);
      survey.set_fagn_pars(lpars);
      LOG_INFO(printf("\nCovariance Matrix:\n"));
      for(int pi=0;pi<pset.covar.size();pi++){
	LOG_INFO(printf("\t["));
	for(int pj=0;pj<pset.covar[pi].size();pj++)
	  LOG_INFO(printf("%6.2lf ",pset.covar[pi][pj]));
        LOG_INFO(printf("]\n"));
      }
      
      LOG_INFO(printf("Final Acceptance Rate: %lf%%\n",metrop.mean_acceptance()*100));
    }
    
    LOG_INFO(printf("\nRunning Best Fit...\n"));
    output=survey.simulate();
    tchi_min = output.chisqr;
    LOG_INFO(printf("Saving Initial Survey\n"));
    saved = survey.save(q.outfile);

    for(pi = 0;pi<q.nparams;pi++)
      means[pi] = pset.best[pi];
    
    for(unsigned long i=1;i<q.nsim;i++){
      rng.gaussian_mv(means,pset.covar,pset.min,pset.max,results);
      for(pi = 0;pi<q.nparams;pi++){
	*(setvars[pi]) = results[pi];
      }
      
      lf.set_params(lpars);
      survey.set_color_exp(CE,ZBC);
      survey.set_fagn_pars(lpars);
      
      output=survey.simulate();
      if(output.chisqr < tchi_min){
	tchi_min = output.chisqr;
	LOG_INFO(printf("Saving new minimum (%lu/%lu)...",i,q.nsim));
	saved = survey.save(q.outfile);
      }
      final_counts.add_link(output.dnds,output.chisqr);
    }
    LOG_INFO(printf("Saved Chi2: %lf (%lf)\n",tchi_min,chi_min));
    
    for(pi=0, p= q.param_inds.begin(); p != q.param_inds.end();pi++,p++)
      parnames[pi] = pnames[*p];
    if(q.vary_cexp)
      parnames[q.cind] = "CEXP";
    if(q.vary_zbc)
      parnames[q.zbcind] = "ZBC";

    LOG_DEBUG(printf("Saving Chains\n"));
    saved &= mcchain.save(q.outfile,parnames.get(),"MCMC Chain Record");
    LOG_DEBUG(printf("Saving Counts\n"));
    saved &= counts.save(q.outfile,countnames,"Simulation Counts");
    LOG_DEBUG(printf("Saving Final Counts\n"));
    saved &= final_counts.save(q.outfile,countnames,"Final Monte Carlo Counts");
  }
  
  saved ? printf("Save Successful") : printf("Save Failed");

  LOG_CRITICAL(printf("\nFitting Complete\n\n"));
  
  return saved ? 0 : 1;
}
