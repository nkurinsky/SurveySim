#include "simulator.h"
#include "functions.h"
#include "constants.h"
#include "cosmo.h"

// find the location of the minimum value
int findMinloc(double vals[], int MAXELS)
{
  int i;
  double min = vals[0];
  int hit=0;

    for (i = 0; i<MAXELS; i++)
    if (min >= vals[i])
      {
	hit=i;
	min = vals[i];
      }
    return hit;
}

sprop::sprop(double z,double f[],double lum, double w){
  redshift = z;
  luminosity = lum;
  weight = w;

  bool do_colors(true);
  for (int i=0;i<3;i++){
    fluxes[i] = f[i];
    if(f[i] < 0.00000000E0)
      do_colors = false;
  }
  if(do_colors){
    c1 = get_color(f[2],f[0]);
    c2 = get_color(f[2],f[1]);
  }
  else{
    c1 = -99.0;
    c2 = -99.0;
  }
}

sprop::sprop(){
  redshift = -1;
  luminosity = -1;
  weight = 0;

  c1 = -99.0;
  c2 = -99.0;
}

void simulator::init_rand(){
  gsl_rng * r;
  const gsl_rng_type * T;
  gsl_rng_default_seed = time(NULL);
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
}

simulator::simulator(double b[],double b_err[],double f_lims[],string obsfile,string sedfile){

  for (int i=0;i<3;i++){
    bands[i] = b[i];
    band_errs[i] = b_err[i];
    flux_limits[i]= f_lims[i];
  }
  
  observations = new obs_lib(obsfile);
  seds = new sed_lib(sedfile);

  //initialize observation portion of histograms
  //possibility here to integrate obs_lib into hist_lib
  diagnostic = new hist_lib();
  double *c1,*c2;
  int osize = observations->get_snum();
  observations->get_all_colors(c1,c2);
  diagnostic->init_obs(c1,c2,osize);
}

void simulator::set_bands(double b[],double b_err[],double f_lims[]){
  for (int i=0;i<3;i++){
    bands[i] = b[i];
    band_errs[i] = b_err[i];
    flux_limits[i]= f_lims[i];
  }
}

void simulator::set_lumfunct(lumfunct *lf){
  if(lf != NULL){
    this->lf = lf;
  }
  else
    cout << "ERROR: NULL Pointer Passed to Simulator" << endl;
}


void simulator::set_sed_lib(string sedfile){
  seds = new sed_lib(sedfile);
}


void simulator::set_obs(string obsfile){
  if(observations != NULL)
    delete observations;
  observations = new obs_lib(obsfile);
  double *c1,*c2;
  int osize = observations->get_snum();
  observations->get_all_colors(c1,c2);
  diagnostic->init_obs(c1,c2,osize);
}

double simulator::simulate(double area, int nz, double dz){
  sources.clear();

  if(seds != NULL){
    //needs to be double checked, should be about done
    static int is,js;
    static double tmpz,vol;
    int lnum = seds->get_lnum();
    double high;
    double scale[nz];
    long nsrcs[nz][lnum];
    double lums[lnum];
    seds->get_lums(lums);

    //setup redshift array    
    double zarray[nz];    
    for (int i=0;i<nz;i++)
      zarray[i]=0.1+i*dz;

    //get the number of sources per L-z bin
    //scale the volume elements by 10^10 to avoid too large a number
    //making a weights array associated with redshift (use volume), aim for 1000 per bin
    //luminosity function dependent scaling (use highest bin, 10 per bin)
    for (is=0;is<nz;is++) {
      tmpz=zarray[is]+dz/2.0;
      vol=(dvdz(tmpz,area)*dz);
      high = lf->get_nsrcs(zarray[is],10.0)*vol;
      scale[is] = 1.0e+4/high;
      for (js=0;js<lnum;js++)
	nsrcs[is][js]= long(scale[is]*vol*lf->get_nsrcs(zarray[is],lums[js]));
    }

    double weights[nz];
    for (int i=0;i<nz;i++)
      weights[i] = 1.e0/scale[i];

    //=========================================================================
    // for each L-z depending on the number of sources, sample the SED and get the appropriate fluxes
    //*************************************************************************
    
    //previous approach was round to nearest template value.
    // We interpolate here, in future will convolve filter
 
    static double b_rest[3],flux_sim[3];
    static double noise[3];
    static sprop *temp_src;
    static int src_iter;
    bool detected = true;

    gsl_rng * r;
    const gsl_rng_type * T;    
    gsl_rng_default_seed = time(NULL);
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    
    //NOTE templates are given in W/Hz
    for (is=0;is<nz;is++){
      b_rest[0]=bands[0]/(1.+zarray[is]);
      b_rest[1]=bands[1]/(1.+zarray[is]);
      b_rest[2]=bands[2]/(1.+zarray[is]);
      printf("\nz = %4.2f, s = %6.2E, s = %6.2E, n: ",zarray[is],scale[is],weights[is]);
      for (js=0;js<lnum;js++){
	printf("{%4.2f, %ld} ",lums[js],nsrcs[is][js]);
	for (src_iter=0;src_iter<nsrcs[is][js];src_iter++){
	  detected = true;
	  for (int i=0;i<3;i++){
	    noise[i]=gsl_ran_gaussian(r,band_errs[i]);
	    //noise[i] = gauss_random(r,nrange,0.0,b_err[i],nsrcs[is][js]); 
	    //this keeps giving me compiler/linker errors
	    flux_sim[i] = seds->get_flux(lums[js],b_rest[i]);
	    flux_sim[i] *= (1/(4.0*M_PI*pow(lumdist(zarray[is])*MPC_TO_METER,2)))/Wm2Hz_TO_mJy;
	    flux_sim[i] += noise[i];
	    if (flux_sim[i] < flux_limits[i]) //reject sources below flux limit
	      detected = false;
	  }	  
	  //check for detectability, if "Yes" add to list
	  if(detected){
	    temp_src = new sprop(zarray[is],flux_sim,lums[js],weights[is]);
	    sources.push_back(*temp_src);
	    delete temp_src;
	  }
	}
      }
    }
    
    //=========================================================================
    // generate diagnostic color-color plots
    //*************************************************************************
    
    double *c1,*c2,*w;
    
    int snum(sources.size());
    c1 = new double[snum];
    c2 = new double[snum];
    w = new double[snum];
    for (int i=0;i<snum;i++){
      c1[i] = sources[i].c1;
      c2[i] = sources[i].c2;
      w[i] = sources[i].weight;
    }
    
    diagnostic->init_model(c1,c2,w,snum);
    chisq=diagnostic->get_chisq();

    delete[] c1;
    delete[] c2;
    delete[] w;

    return chisq;
  }
  else{
    cout << "ERROR: NULL Model Library" << endl;
    return -1;
  }
}


bool simulator::save(string outfile){
  return diagnostic->write_fits(outfile);
}

simulator::~simulator(){
  if(observations != NULL)
    delete observations;
  if(seds != NULL)
    delete seds;  
  if(diagnostic != NULL)
    delete diagnostic;
}

void simulator::reset(){
  sources.clear();
}