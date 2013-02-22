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

sprop::sprop(int m,double z,double b[],double f[],double lum){
  mod_num = m;
  redshift = z;
  luminosity = lum;
  bool do_colors(true);
  for (int i=0;i<3;i++){
    bands[i] = b[i];
    fluxes[i] = f[i];
    if((b[i] < 0.0000000E0) or (f[i] < 0.00000000E0))
      do_colors = false;
  }
  if(do_colors){
    c1 = log10(f[2]/f[0]); // /log10(b[2]/b[0]);
    c2 = log10(f[2]/f[1]); // /log10(b[2]/b[1]);
  }
  else{
    c1 = -99.0;
    c2 = -99.0;
  }
}

sprop::sprop(){
  mod_num = -1;
  redshift = -1;
  luminosity = -1;
  /*
  for (int i=0;i<3;i++){
    fluxes[i] = -1;
    bands[i] = -1;
  }
  */
  c1 = -99.0;
  c2 = -99.0;
}

/*
simulator::simulator(){
  //distributions = NULL;
  luminosity_function = NULL;
  models = NULL;
  //  observations = NULL;
  diagnostic = NULL;
  //distribution_size = 1000;

  redshift_range[0] = 0.0;
  redshift_range[1] = 10.0;

  for (int i=0;i<3;i++){
    bands[i] = -1;
    band_errs[i] = -1;
    flux_limits[i]= -1;
  }
}
*/

void simulator::init_rand(){
  gsl_rng * r;
  const gsl_rng_type * T;
  gsl_rng_default_seed = time(NULL);
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
}

simulator::simulator(double b[],double b_err[],double area,double f_lims[],double lpars[],string modfile,string obsfile,string sedfile){

  cout<<"from inside simulator"<<endl;

  //========================================================================
  //1st step: define the lambda and luminosity array 
  //these are really set up when the sf_templates.fits is set
  //*************************************************************************

  FITS *pInfile;
  
  pInfile = new FITS(sedfile,Read);

  std::valarray<double> seds;
  PHDU& image = pInfile->pHDU();
  image.read(seds); //store contents of image into the array "seds"

  int nlambda = image.axis(0);
  int nlum = image.axis(1)-1; //since the first row here is the lambda array

  int fi;
  double lambda[nlambda];
  //for (int fi=0;fi<nlambda;fi++){
  //  lambda[fi]=seds[fi];
  //}

 //read in the total IR luminosities associated with the different z=0 SED templates
  double lumarray[nlum];
  double temp;
  HDU& header = pInfile->pHDU();
  
  for (int i=0;i<nlum;i++) {
    std::string a = "ROW";
    a+= static_cast<ostringstream*>( &(ostringstream() << i+1) )->str();
    header.readKey(a,temp);
    lumarray[i]=temp;  
  }
 
  //=========================================================================
  //2nd step: setup redshift array
  //*************************************************************************

  int nz=10;
  double dz=0.5;
  double zarray[nz];

  for (int fi=0;fi<nz;fi++)
    zarray[fi]=0.1+fi*dz;

  double seds0[nlambda][nlum];

  for (int fi=0;fi<nlambda;fi++){
    lambda[fi]=seds[fi];
    for (int fj=0;fj<nlum;fj++){
      seds0[fi][fj]=seds[fi+nlambda*(fj+1)];
    }
  }
  
  //=========================================================================
  //3rd step: determine the # of sources in each L-z bin, depending on the LF model
  //*************************************************************************

  lumfunct lf; //
  
  lf.set_phi0 (lpars[0]);
  lf.set_L0 (lpars[1]);
  lf.set_alpha (lpars[2]);
  lf.set_beta (lpars[3]);
  lf.set_p (lpars[4]);
  lf.set_q (lpars[5]);

  int is,js;
  double tmpz,vol;
  double nsrcs[nz][nlum];

  //get the number of sources per L-z bin
  //scale the volume elements by 10^5 to avoid too large a number
  for (is=0;is<nz;is++) 
    {
      tmpz=(zarray[is+1]+zarray[is])/2.0;
      vol=(dvdz(tmpz,area)*dz)/(1.e+05);
      for (js=0;js<nlum;js++)
	{
	  nsrcs[is][js]=vol*lf.get_nsrcs(zarray[is],lumarray[js]);
	  //cout<<nsrcs[is][js]<<endl;
	}
    }
  
  //=========================================================================
  //4th step: for each L-z depending on the number of sources, sample the SED and get the appropriate fluxes
  //*************************************************************************

  cout<<"Bands of interest"<<endl;
  cout<<b[0]<<endl;
  cout<<b[1]<<endl;
  cout<<b[2]<<endl;

  int hit1,hit2,hit3;
  double testlam[nlambda];
  double b_rest[3],flux_sim[3];

  //note this is not very precise (see the actual b-rest and lambda[hit] values)
  //ultimately want to convolve with filter but this is quick for now. 

  static double nrange[2];
  //double ** noise = new double*[3];
  double noise[3];

  sources.clear();

  static sprop *temp_src;

  int src_iter;
  
  gsl_rng * r;

  const gsl_rng_type * T;
  //gsl_rng_env_setup();

  gsl_rng_default_seed = time(NULL);
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  //NOTE templates are given in W/Hz
  for (is=0;is<nz;is++)
    {
      //cout<<zarray[is]<<endl;
      b_rest[0]=b[0]/(1.+zarray[is]);
      for (fi=0;fi<nlambda;fi++) testlam[fi]=abs(lambda[fi]-b_rest[0]);
      hit1=findMinloc(testlam,nlambda);
      b_rest[1]=b[1]/(1.+zarray[is]);
      for (fi=0;fi<nlambda;fi++) testlam[fi]=abs(lambda[fi]-b_rest[1]);
      hit2=findMinloc(testlam,nlambda);
      b_rest[2]=b[2]/(1.+zarray[is]);
      for (fi=0;fi<nlambda;fi++) testlam[fi]=abs(lambda[fi]-b_rest[2]);
      hit3=findMinloc(testlam,nlambda);
      for (js=0;js<nlum;js++)
	{
	  //printf("L-z bin: %lf %lf \n",zarray[is],lumarray[js]);
	  //cout<<"Number of sources"<<endl;
	  //cout<<nsrcs[is][js]<<endl;
	  
	  
	  for (src_iter=0;src_iter<nsrcs[is][js];src_iter++)
	    {
	      for (int i=0;i<3;i++){
  //	    nrange[0] = b_err[i]*-5.0;
  //    nrange[1] = b_err[i]*5.0;
	    //noise[i]=0;
	    noise[i]=gsl_ran_gaussian(r,b_err[i]);
	    //noise[i] = gauss_random(r,nrange,0.0,b_err[i],nsrcs[is][js]); //this keeps giving me compiler/linker errors
	      }

	      //cout<<src_iter<<endl;
	      
	      flux_sim[0]=(seds0[hit1][js]/(4.0*M_PI*pow(lumdist(zarray[is])*MPC_TO_METER,2)))/Wm2Hz_TO_mJy;
	      flux_sim[0]+=noise[0]; //[src_iter];
	      flux_sim[1]=(seds0[hit2][js]/(4.0*M_PI*pow(lumdist(zarray[is])*MPC_TO_METER,2)))/Wm2Hz_TO_mJy;
	      flux_sim[1]+=noise[1]; //[src_iter];
	      flux_sim[2]=(seds0[hit3][js]/(4.0*M_PI*pow(lumdist(zarray[is])*MPC_TO_METER,2)))/Wm2Hz_TO_mJy;
	      flux_sim[2]+=noise[2]; //[src_iter];
	      //check for detectability, if "Yes" add to list
	      //	      if(flux_sim[0] >= 8.0) //include real flux limits but had errors for some reason
              cout<<f_lims[0]<<endl;
              if(flux_sim[0] >=f_lims[0])
	      {
		//cout<<flux_sim[0]<<endl;
		temp_src = new sprop(0.0,zarray[is],b,flux_sim,lumarray[js]);
		//cout<<temp_src->redshift<<endl;
		//print,flux_sim[0]
		sources.push_back(*temp_src);
		//cout<<sources.redshift[src_iter]<<endl;
		delete temp_src;
	      }
	    }
	}
    }

  //cout<<"Size of sources vector: "<<sources.redshift()<<endl;

 //=========================================================================
 //5th step: generate diagnostic color-color plots
 //*************************************************************************

  double *mc1,*mc2;

  int snum(sources.size());
  mc1 = new double[snum];
  mc2 = new double[snum];
  for (int i=0;i<snum;i++){
    mc1[i] = sources[i].c1;
    mc2[i] = sources[i].c2;
  }

  //the observed and model colors
  //double *c1,*c2,*mc1,*mc2;
  

  //int osize = observations->get_snum();
  int msize = nz*nlum; //distributions->size();

  //cout<<osize<<endl;
  cout<<" The n_z x n_lum size: "<<endl;
  cout<<msize<<endl;

    if (diagnostic == NULL)
      diagnostic = new hist_lib;

    //for(src_iter=0;src_iter<nsrcs[0][0];src_iter++)
    //{
    //cout<<sources.redshift[src_iter]<<endl;
    //}
    //cout<<sources.fluxes[0]<<endl;
    //cout<<sources.c1[0]<<endl;
    //mc1=sources.c1;
    //mc2=sources.c2;
  
  //    observations->get_all_colors(c1,c2);
    models->get_all_colors(mc1,mc2);

    // struct sprop model;
    //cout<<sources.redshift[0]<<endl;
    //diagnostic->init_obs(c1,c2,osize);
    diagnostic->init_model(mc1,mc2,msize);
    m_chi2=diagnostic->get_chisq();

    //delete[] c1;
    //delete[] c2;
    delete[] mc1;
    delete[] mc2;
 
    
    // return diagnostic->get_chisq();

  //tmpz=zarray[0];
  //tmpl=lumarray[0];

  // sed_models model_lib(sedfile);

  //fluxes=get_fluxes(seds,b,tmpz,tmpl);
  
  /*
  //distributions = NULL;
  diagnostic = NULL;
  //distribution_size = 1000;

  redshift_range[0] = 0.0;
  redshift_range[1] = 10.0;

  if(lf != NULL)
    luminosity_function = lf;
  else
    cout << "ERROR: NULL Pointer Passed to Simulator" << endl;
  
  models = new model_lib(modfile);
  observations = new obs_lib(obsfile);
  
  for (int i=0;i<3;i++){
    bands[i] = b[i];
    band_errs[i] = b_err[i];
    flux_limits[i]= f_lims[i];
  }
  */
}

void simulator::set_bands(double b[],double b_err[],double f_lims[]){
  for (int i=0;i<3;i++){
    bands[i] = b[i];
    band_errs[i] = b_err[i];
    flux_limits[i]= f_lims[i];
  }
}




/*
void simulator::set_lumfunct(lumfunct *lf){
  if(lf != NULL)
    luminosity_function = lf;
  else
    cout << "ERROR: NULL Pointer Passed to Simulator" << endl;
}
*/
//void simulator::set_model_lib(string modfile){
//  models = new model_lib(modfile);
//}

/*
void simulator::set_obs_lib(string obsfile){
  observations = new obs_lib(obsfile);
}
*/

//void simulator::set_zrange(double zrange[]){
// if(zrange[0] > 0.0)
//    redshift_range[0] = zrange[0];
//  if(zrange[1] < 10.0)
//    redshift_range[1] = zrange[1];
//}

/*
double simulator::simulate(const gsl_vector *v,void *params){
  reset();

  if(models != NULL){
    //double mean[3];
    //double sigma[3];
    
    //fixed_params *p = (fixed_params *)params;
    
    //mean[0] = gsl_vector_get(v,0);
    //sigma[0] = abs(gsl_vector_get(v,1));
    
    //for (int i=1;i<=p->pnum;i++){
    //  mean[i] = p->p[i-1].mean;
    //  sigma[i] = p->p[i-1].sigma;
    // }
    
    //distributions = new sim(models);
    //distributions->initialize_dists(mean,sigma,distribution_size,redshift_range);

    if(bands[0] != -1){
      //      distributions->simulate(bands,band_errs,flux_limits,*luminosity_function);
    }
    else{
      cout << "ERROR: Bands Uninitialized" << endl;
      return -1;
    }
    
    double *c1,*c2,*mc1,*mc2;
    int osize = observations->get_snum();
    //int msize = distributions->size();
    diagnostic = new hist_lib;

    observations->get_all_colors(c1,c2);
    //distributions->get_all_colors(mc1,mc2);

    diagnostic->init_obs(c1,c2,osize);
    int msize=1; //PLACEHOLDER!!!!
    diagnostic->init_model(mc1,mc2,msize);
    
    return diagnostic->get_chisq();
  }
  else{
    cout << "ERROR: NULL Model Library" << endl;
    return -1;
  }
}
*/

bool simulator::save(string outfile){
  return diagnostic->write_fits(outfile);
}

simulator::~simulator(){

  //if(distributions != NULL)
  //  delete distributions;

  //  if(models != NULL)
  //  delete models;
  
  //  if(observations != NULL)
  //delete observations;

  if(diagnostic != NULL)
    delete diagnostic;
}

void simulator::reset(){
  
  if(diagnostic != NULL)
    delete diagnostic;

  //if(distributions != NULL)
  //  delete distributions;
}




