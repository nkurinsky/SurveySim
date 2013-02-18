
#include "functions.h"

using namespace std;

gsl_rng * r;  /* global generator */

int main(){

const gsl_rng_type * T;

 gsl_rng_default_seed=time(NULL);
    T=gsl_rng_default;
    r=gsl_rng_alloc(T);

  double testg;

  static double nr[2];
  nr[0]=-1;
  nr[1]=1;

  int gsize=1; //the size for gauss_random

  double pmean,qmean;
  pmean=0.0;
  qmean=0.0;

  int npar=2;
  double dp[npar];
  dp[0]=0.5;
  dp[1]=1.0;

  //printf ("generator type: %s\n", gsl_rng_name (r));
  //    printf ("seed = %lu\n", gsl_rng_default_seed);
  //     printf ("first value = %lu\n", gsl_rng_get (r));

  //double gsl_ran_gaussian (const gsl_rng * r, double sigma)
  testg=gsl_ran_gaussian(r,dp[0]);
  //the gauss_random keeps giving me errors
  //  testg=gauss_random(r,nr,pmean,dp[0],gsize);
  cout<<testg<<endl;

  gsl_rng_free (r);
  return 0;
}
