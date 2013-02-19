//****************************************************************************
// written by Anna Sajina 09/28/12
// The purpose of this program is to read-in the parameters set in the idl wrapper (simulation.pro) and to use MCMC to determine the best-fit set of parameters. // Specifically it runs through a number of iterations (runs) in each stage calling on simulator.cpp to compute the chi2 for a given set of parameters
// then the metrop algorithm is invoked to decide whether or not a given trial run is accepted. If so it is added to a chain of mcmc values.
//****************************************************************************

#include "simulator.h"
#include "functions.h"

using namespace std;

gsl_rng * r;  /* global generator */

//here de is effectively Delta(Energy), and t is the temperature
//as in the boltzman factor exp(-de/t)
bool metrop(double de,double t){
  
  bool ans;
  
  const gsl_rng_type * T;
  gsl_rng_env_setup();   

  gsl_rng_default_seed = time(NULL);
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  //want uniform random number from 0-1
  double itemp=gsl_rng_uniform(r);
  //cout<<itemp;

  if ((de < 0.0) or (itemp < exp(-de/t))){
      ans=true;
    } else ans=false;
  
  return (ans);
}

//note here we have arguments that are being passed to this (as in by the idl widget), but never seem to do anything with them
int main(int argc,char** argv){
  
  long int runs=30; //temporary, of course in the end this will be much longer perhaps 100,000-500,000 
  
  const gsl_rng_type * T;
  
  int npar=2; 
  double p_o[npar],dp[npar],p_min[npar],p_max[npar]; //the initial guesses of the parameters, the width of the proposal distribution and the acceptable min/max range
  double chain[npar][runs];

  //Note, these have to come from the widget, not be hardwired here!
  string outfile("/Users/annie/students/noah_kurinsky/Fitting/output.fits");
  string modfile("/Users/annie/students/noah_kurinsky/Fitting/model.fits");
  string sedfile("/Users/annie/students/noah_kurinsky/Fitting/sf_templates.fits");
  string obsfile("/Users/annie/students/noah_kurinsky/Fitting/observation.fits");
  
  if(argc > 1)
    outfile = argv[1];
  
  int pnum; 
  
  FITS *pInfile,*pInfile2;
  
  pInfile = new FITS(modfile,Read);
  pInfile2 = new FITS(obsfile,Read);
  
  //reading primary header from fits file
  HDU& params_table = pInfile->pHDU();
  HDU& obs_table = pInfile2->pHDU();

  //get number of parameters and initialize dynamic memory
  params_table.readKey("P_NUM",pnum);

  string pi[] = {"0","1","2","3","4","5"};

  int bands=3;
  double bs[bands],errs[bands],flims[bands];

  for(int i=0;i<3;i++){
    obs_table.readKey("WAVE_"+pi[i+1],bs[i]);
    obs_table.readKey("W"+pi[i+1]+"_FMIN",flims[i]);
    obs_table.readKey("W"+pi[i+1]+"_FERR",errs[i]);
    //bs[i]*=1.e+06; //convert to microns -- Not needed, provide in microns to begin with!
    //flims[i]*=1e-3; //this is to turn it in Janskys, but I think lets work in mJy anyway
    //errs[i]*=1e-3;
  }
  
//=================================================================
// Before doing any fits, initialize the observed diagnostic histogram
// best to do it here, so don't have to re-do it everytime we call on
// simulator
//------------------------------------------------------------------

  //double *c1,*c2;
  //  hist_lib hist_lib(obsfile,bands);
  //obs_lib obs;
  //observations.get_all_colors(c1,c2);
  //int osize;
  //obs_lib get_snum();
  //hist_lib();
  //hist_lib(c1,c2,osize);
    //hist_lib


//=================================================================  
//Read-in Luminosity Function Parameters
//-----------------------------------------------------------------

//this array holds phi0,L0,alpha,beta,p, and q
  double lpars[6];
  double temp;
  params_table.readKey("PHI0",temp);
  lpars[0]=temp;
  params_table.readKey("L0",temp);
  lpars[1]=temp;
  params_table.readKey("ALPHA",temp);
  lpars[2]=temp;
  params_table.readKey("BETA",temp);
  lpars[3]=temp;
  params_table.readKey("P",temp);
  lpars[4]=temp;
  params_table.readKey("Q",temp);
  lpars[5]=temp;
  
  //note that the dp values here are the widths of the "proposal distribution"
  //the smaller they are the slower will converge onto the right answer
  //the bigger they are, the less well sampled the probability distribution will be and hence the less accurate the final answer
  //needs to find the "goldilocks" zone here but that's very problem specific so requires trial and error as we actually start fitting data.
  p_o[0]=lpars[4]; 
  dp[0]=0.3;
  p_min[0]=0.0;
  p_max[0]=7.0;
  
  p_o[1]=lpars[5];
  dp[1]=0.3;
  p_min[1]=0.0;
  p_max[1]=5.0;
  
  printf("Initial p: %5.3f, and q: %5.3f\n",p_o[0],p_o[1]);
  fixed_params ps;
  ps.pnum = pnum;

  delete pInfile;
  delete pInfile2;
  
  //loop trough the mcmc runs
  // at each step, from some initial paramer value "p" call on the simulator.cpp //program to evaluate the chi2 for that particular color-color plot. Then use the metrop algorithm (below) to decide whether or not to keep a particular guess
  //the chain results (including all accepted guesses) will provide us with the posterior probability distributions for the fitted parameters.

  double chi_min=1000.0; //the initial chi2_min value, this is iterated each time a better value is found
  double previous=chi_min;
  double trial;
  int minlink;
  long acceptot=0;
  bool ans;

  double p0_rng[runs];
  double q0_rng[runs];

  gsl_rng_default_seed=time(NULL);
  T=gsl_rng_default;
  r=gsl_rng_alloc(T);
  
  double area;
  double test_chi2;

  //note need to be able to pass the survey area down from the widget!!!
  //this is necessary for correct cosmological volume determinations
  //hence correct number count predictions
  //this value of area here is just a placeholder
  area=pow((1.4*M_PI/180.0),2);

  double t=100; // the temperature, keep fixed for now, but can also try annealing in the future, this has similar effect as the width of the proposal distribution, as determines how likely a far flung guess is of being accepted (see metrop function on top)

    //REMOVE THIS, TESTING PURPOSES ONLY as the value in params.save is currently wrong
    p_o[0]=6.0;

  for (int i=0;i<runs;i++){
    p0_rng[i]=gsl_ran_gaussian(r,dp[0]);
    q0_rng[i]=gsl_ran_gaussian(r,dp[1]);
    temp=p_o[0]+p0_rng[i];  
    if((temp >= p_min[0]) && (temp <= p_max[0])) lpars[4]=temp;
    temp=p_o[1]+q0_rng[i];
    if((temp >= p_min[1]) && (temp <= p_max[1])) lpars[5]=temp;
    // this should be commented out when the runs number gets longer as it will slow things down alot to have to display each guess
    cout<<i+1<<" "<<lpars[4]<<" "<<lpars[5]<<endl; //check to see if sensible guesses, need to also do some test the randomness at some point

    /*
    simulator tester(bs,errs,area,flims,lpars,modfile,obsfile,sedfile);
    test_chi2=tester.model_chi2();

    cout<<"Model chi2 :"<<test_chi2<<endl;

    trial=test_chi2;
    double de=trial-chi_min;
    if(trial < chi_min){
      chi_min=trial;
      minlink=i;
    }
    ans=metrop(de,t);
    if(ans=true){
      acceptot++;
      previous=trial;
      //update mcmc chain with accepted values
	chain[0][i]=lpars[4];
	chain[1][i]=lpars[5];
    }
    */
  }

  gsl_rng_free (r);

  return 0;
}


