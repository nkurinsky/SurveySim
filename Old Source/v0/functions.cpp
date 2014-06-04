#include "functions.h"

//Generates uniform distribution using GNU Standard algorithm (mt19937)
double * random(gsl_rng * r,double range[],int size){
  double * temp = new double[size];
  static double itemp;
  for (int i=0;i<size;i++){ //generates array of specific size
    itemp = gsl_rng_uniform(r); //GNU function (ch 18)
    itemp *= (range[1]-range[0]);
    itemp += range[0];
    temp[i] = itemp;
  }
  return temp;
}

//Computes random gaussian distribution with given properties
double * gauss_random(gsl_rng * r,double range[],double mean,double sigma,int size){
  double * temp = new double[size];
  static double itemp;
  for (int i=0;i<size;i++){
    do{
      itemp = gsl_ran_gaussian(r,sigma);
      itemp += mean;
    }while((itemp > range[1]) or (itemp < range[0]));
    temp[i] = itemp;
  }
  return temp;
}

//Computes spectral color in log space from flux and band inputs
double get_color(double f1,double f2){
  return log10(f1/f2);
}
