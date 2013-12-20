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

//Array size definitions (most likely will never change)
#define BANDS 3
//For values to be passed in by widget
#define NPAR 2
#define NCHAIN 5
#define TMAX 20.00
#define IDEALPCT 0.25
#define BURN_STEP 10
#define CONV_STEP 20

using namespace std; 

gsl_rng * r;  // global generator

int main(int argc,char** argv){
  
  if(argc < 4){
    printf("%s","ERROR: Invalid number of arguments.\n");
    printf("%s","Calling sequence should be \"fit obsfile modfile sedfile [output]\"\n");
    return 1;}

  printf("\nMCMC Fitting Code Compiled: %s\n",__DATE__);
  
  //File names passed in by Widget
  string outfile("output.fits");
  string obsfile(argv[1]);
  string modfile(argv[2]);
  string sedfile(argv[3]);
  //If outfile specified as argument, change from default
  if(argc > 4)
    outfile = argv[4];

  unsigned long runs; 
  int nz,ns;
  double area, dz, zmax, zmin, rtemp;
  products output;
  double bs[BANDS],errs[BANDS],flims[BANDS];
  string pi[] = {"0","1","2","3","4","5"};
  //these arrays holds phi0,L0,alpha,beta,p, and q as well as range and fixed
  double lpars[6],lfix[6],lmax[6],lmin[6],ldp[6],cexp[5];
  // the initial guesses of the parameters, the width of the proposal distribution 
  // and the acceptable min/max range
  double p_o[NPAR],dp[NPAR],p_min[NPAR],p_max[NPAR]; 

  const gsl_rng_type * T;

  FITS *pInfile,*pInfile2;
  pInfile = new FITS(modfile,Read);
  pInfile2 = new FITS(obsfile,Read);
  
  //reading primary header from fits file
  HDU& params_table = pInfile->pHDU();
  HDU& obs_table = pInfile2->pHDU();

  params_table.readKey("RUNS",rtemp);
  runs = (unsigned long) rtemp;
  params_table.readKey("ZMIN",zmin);
  params_table.readKey("ZMAX",zmax);
  params_table.readKey("DZ",dz);
  params_table.readKey("AREA",area);

  //this is necessary for correct cosmological volume determinations
  //hence correct number count predictions
  area*=pow((M_PI/180.0),2.0); //in steradians from sq.deg.

  //compute redshift bin number from FITS values
  nz = (zmax-zmin)/dz;
  printf("\nNR: %lu, NZ: %i, DZ: %f\n",runs,nz,dz);

  printf("%s\n","Bands:");
  for(int i=0;i<BANDS;i++){
    obs_table.readKey("WAVE_"+pi[i+1],bs[i]); //should already be in microns
    obs_table.readKey("W"+pi[i+1]+"_FMIN",flims[i]); //should be in mJy
    obs_table.readKey("W"+pi[i+1]+"_FERR",errs[i]);
    printf("%fl\t%fl\t%fl\n",bs[i],flims[i],errs[i]);
  }
  
  //=================================================================  
  //Read-in Luminosity Function Parameters
  //-----------------------------------------------------------------

  //read parameter values
  params_table.readKey("PHI0",lpars[0]);
  params_table.readKey("L0",lpars[1]);
  params_table.readKey("ALPHA",lpars[2]);
  params_table.readKey("BETA",lpars[3]);
  params_table.readKey("P",lpars[4]);
  params_table.readKey("Q",lpars[5]);
  //read parameter value fitting boolean
  params_table.readKey("PHI0_FIX",lfix[0]);
  params_table.readKey("L0_FIX",lfix[1]);
  params_table.readKey("ALPHA_FI",lfix[2]);
  params_table.readKey("BETA_FIX",lfix[3]);
  params_table.readKey("P_FIX",lfix[4]);
  params_table.readKey("Q_FIX",lfix[5]);
  //read parameter minimum values
  params_table.readKey("PHI0_MIN",lmin[0]);
  params_table.readKey("L0_MIN",lmin[1]);
  params_table.readKey("ALPHA_MI",lmin[2]);
  params_table.readKey("BETA_MIN",lmin[3]);
  params_table.readKey("P_MIN",lmin[4]);
  params_table.readKey("Q_MIN",lmin[5]);
  //read parameter maximum values
  params_table.readKey("PHI0_MAX",lmax[0]);
  params_table.readKey("L0_MAX",lmax[1]);
  params_table.readKey("ALPHA_MA",lmax[2]);
  params_table.readKey("BETA_MAX",lmax[3]);
  params_table.readKey("P_MAX",lmax[4]);
  params_table.readKey("Q_MAX",lmax[5]);
  //read parameter sigma values
  params_table.readKey("PHI0_DP",ldp[0]);
  params_table.readKey("L0_DP",ldp[1]);
  params_table.readKey("ALPHA_DP",ldp[2]);
  params_table.readKey("BETA_DP",ldp[3]);
  params_table.readKey("P_DP",ldp[4]);
  params_table.readKey("Q_DP",ldp[5]);

  //read color_exp values
  params_table.readKey("CEXP",cexp[0]);
  params_table.readKey("CEXP_FIX",cexp[1]);
  params_table.readKey("CEXP_MIN",cexp[2]);
  params_table.readKey("CEXP_MAX",cexp[3]);
  params_table.readKey("CEXP_DP",cexp[4]);
  
  /*
    note that the dp values here are the widths of the "proposal distribution"
    the smaller they are the slower will converge onto the right answer
    the bigger they are, the less well sampled the probability distribution 
    will be and hence the less accurate the final answer
    needs to find the "goldilocks" zone here but that's very problem specific 
    so requires trial and error as we actually start fitting data. 
  */
  p_o[0] = lpars[4]; 
  p_o[1] = lpars[5];
  dp[0] = ldp[4];
  dp[1] = ldp[5];
  p_min[0] = lmin[4];
  p_min[1] = lmin[5];
  p_max[0] = lmax[4];
  p_max[1] = lmax[5];

  lumfunct lf;
  lf.set_phi0(lpars[0]);
  lf.set_L0(lpars[1]);
  lf.set_alpha(lpars[2]);
  lf.set_beta(lpars[3]);
  lf.set_p(lpars[4]);
  lf.set_q(lpars[5]);
  
  printf("Initial p: %5.3f, and q: %5.3f\n",p_o[0],p_o[1]);
  printf("p range: %5.3f - %5.3f, sigma: %5.3f\n",p_min[0],p_max[0],dp[0]);
  printf("q range: %5.3f - %5.3f, sigma: %5.3f\n",p_min[1],p_max[1],dp[1]);
  printf("Color Evolution: %5.3f",cexp[0]);

  delete pInfile;
  delete pInfile2;

  /*
    Loop trough the mcmc runs
    At each step, from some initial paramer value "p" call on the simulator.cpp 
    program to evaluate the chi2 for that particular color-color plot. 
    Then use the metrop algorithm (below) to decide whether or not to keep a particular guess
    The chain results (including all accepted guesses) will provide us with the posterior 
    probability distributions for the fitted parameters.
  */

  //initialize simulator
  simulator survey(bs,errs,flims,obsfile,sedfile); 
  survey.set_color_exp(cexp[0]); //color evolution
  survey.set_lumfunct(&lf);

  //the initial chi2_min value
  //this is iterated each time a better value is found
  double chi_min=1.0E+4; 
  double trial;
  double prng[NPAR][runs];
  unsigned long i,m,p;

  // the temperature, keep fixed for now, but can also try annealing 
  // in the future, this has similar effect as the width of the proposal
  // distribution, as determines how likely a far flung guess is of being 
  // accepted (see metrop function)
  //double tmc=10.00; //to distinguish it from the random T
  double temp;
  double ptemp[NCHAIN][NPAR];
  double pbest[NPAR];
  double pcurrent[NCHAIN][NPAR];

  gsl_rng_default_seed=time(NULL);
  T=gsl_rng_default;
  r=gsl_rng_alloc(T);

  //the flux array is logarithmic in steps of 0.3dex 
  //for now lets just use one band (here 250um) although of course might be nice to keep the rest at some point, but one step at a time.
  ns=8;
  survey.set_size(area,dz,zmin,nz,ns);
  MCChains mcchain(NCHAIN,NPAR,runs);
  MetropSampler metrop(NCHAIN,TMAX,IDEALPCT,r);

  for(m=0;m<NCHAIN;m++){    
    ptemp[m][0] = pcurrent[m][0] = lpars[4];
    ptemp[m][1] = pcurrent[m][1] = lpars[5];
    switch(m){
    case 1:
      ptemp[1][0] += (p_max[0]-ptemp[1][0])/2.0;
      pcurrent[1][0] = ptemp[1][0];
      break;
    case 2:
      ptemp[2][0] += (p_min[0]-ptemp[2][0])/2.0;
      pcurrent[2][0] = ptemp[2][0];
      break;
    case 3:
      ptemp[3][1] += (p_max[1]-ptemp[3][1])/2.0;
      pcurrent[3][1] = ptemp[3][1];
      break;
    case 4:
      ptemp[4][1] += (p_min[1]-ptemp[4][1])/2.0;
      pcurrent[4][1] = ptemp[4][1];
      break;
    default:
      break;
    }
  }
  
  printf("\n ---------- Beginning MC Burn-in Phase ---------- \n");
  for (i=0;i<(runs/10);i++){
    for (m=0;m<NCHAIN;m++){
      for(p=0;p<NPAR;p++){
	prng[p][i]=gsl_ran_gaussian(r,dp[p]);
	temp=pcurrent[m][p]+prng[p][i];
	if((temp >= p_min[p]) && (temp <= p_max[p])) ptemp[m][p]=temp;
      }
      
      printf("%lu %lu - %lf %lf : ",(i+1),(m+1),ptemp[m][0],ptemp[m][1]);
      
      //this not yet generalized
      lf.set_p(ptemp[m][0]);
      lf.set_q(ptemp[m][1]);
      output=survey.simulate();
      trial=output.chisqr;
      printf("Iteration Chi-Square: %lf\n",trial);
      
      if(trial < chi_min){
	printf(" -- Current Best Trial --");
	chi_min=trial;
	for(p=0;p<NPAR;p++)
	  pbest[p]=ptemp[m][p];
      }

      if(metrop.accept(m,trial)){
	for(p=0;p<NPAR;p++)
	  pcurrent[m][p]=ptemp[m][p];
      }      
    }
    if(((i+1) % BURN_STEP) == 0)
      if(not metrop.anneal())
	i = runs;
  }
  
  for(m=0;m<NCHAIN;m++)
    for(p=0;p<NPAR;p++)
      ptemp[m][p] = pcurrent[m][p] = pbest[p];
  
  metrop.reset();
  
  printf("\n\n ---------------- Fitting Start ---------------- \n Total Run Number: %ld\n\n",runs);
  for (i=0;i<runs;i++){
    for (m=0;m<NCHAIN;m++){
      for(p=0;p<NPAR;p++){
	prng[p][i]=gsl_ran_gaussian(r,dp[p]);
	temp=pcurrent[m][p]+prng[p][i];
	if((temp >= p_min[p]) && (temp <= p_max[p])) ptemp[m][p]=temp;
      }
      
      printf("%lu %lu - %lf %lf : ",(i+1),(m+1),ptemp[m][0],ptemp[m][1]);
      lf.set_p(ptemp[m][0]);
      lf.set_q(ptemp[m][1]);
      output=survey.simulate();
      trial=output.chisqr;
      printf("Model chi2: %lf",trial);
      
      mcchain.add_link(m,ptemp[m],trial);
      
      if(metrop.accept(m,trial)){
	printf(" -- Accepted\n");
	for(p=0;p<NPAR;p++)
	  pcurrent[m][p]=ptemp[m][p];
      }
      else
	printf(" -- Rejected\n");
    }
    
    if(((i+1) % CONV_STEP) == 0){
      printf("Checking Convergence\n");
      metrop.anneal();
      if(mcchain.converged()){
	printf("Chains Converged!\n");
	i = runs;}
      else
	printf("Chains Have Not Converged\n");
    } 
  }
  
  mcchain.get_best_link(pbest,chi_min);
  lf.set_p(pbest[0]);
  lf.set_q(pbest[1]);
  printf("\nRe-Running Best Fit...\n");
  output=survey.simulate();
  printf("Model chi2: %lf\n",output.chisqr);
  printf("Acceptance Rate: %lf%%\n",metrop.mean_acceptance());
  
  string parnames[] = {"p0","q0"};
  bool saved;
  saved = survey.save(outfile);
  saved &= mcchain.save(outfile,parnames);
  if(saved)
    printf("Save Successful\n");
  else
    printf("Save Failed\n");
  
  gsl_rng_free (r);
  
  printf("\nReturning from fitting routine\n\n");
  
  return 0;
}
