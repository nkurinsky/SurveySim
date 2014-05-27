#include "obs_lib.h"

obs::obs(){
  c1 = -99;
  c2 = -99;
  interp_init = false;
}

obs::obs(double * f,double *bands){
  bool do_colors = true;
  for (int i=0;i<3;i++){
    fluxes[i] = f[i];
    if((bands[i] <= 0.0000000E0) or (f[i] <= 0.00000000E0))
      do_colors = false;
  }

  if(do_colors){
    c1 = get_color(fluxes[2],fluxes[0]);
    c2 = get_color(fluxes[2],fluxes[1]);
  }
  else{
    c1 = -99.0;
    c2 = -99.0;
  }

  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline,3);
  gsl_spline_init(spline,bands,f,3);
  interp_init = true;
}

double obs::get_flux(double band){
  return gsl_spline_eval(spline,band,acc);
}

void obs::get_colors(double &c1,double &c2){
  c1 = this->c1;
  c2 = this->c2;
}

obs::~obs(){
  if(interp_init){
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
}

obs_lib::obs_lib(string fitsfile){
  //initialize FITS input
  FITS* pInfile;
  pInfile = new FITS(fitsfile,Read);
  
  //reading primary header from fits file
  HDU& table = pInfile->pHDU();
  //extracting conents of image from fits file
  std::valarray<double> contents;
  PHDU& image = pInfile->pHDU();
  image.read(contents);

  static string ti[3] = {"1","2","3"};
  double btemp,fetemp,fltemp;
  observations.resize(image.axis(1));

  if(image.axis(0) != 3){
    printf("%s \n","ERROR: (Obs_Lib) Incorrect Number of Bands");
  }
  else{

    for (int i=0;i<3;i++){
      table.readKey("WAVE_"+ti[i],btemp);
      table.readKey("W"+ti[i]+"_FERR",fetemp);
      table.readKey("W"+ti[i]+"_FMIN",fltemp);
      bands[i]=btemp;
      ferr[i]=fetemp;
      flim[i]=fltemp;
    }
    
    obs *new_obs;
    double fluxes[3];
    
    for (unsigned int i=0;i<observations.size();i++){
      for (int j=0;j<3;j++)
	fluxes[j]=contents[i*3+j];
      new_obs = new obs(fluxes,bands);
      observations[i] = new_obs;
    }
  }
  
  delete pInfile;
}
  
double obs_lib::get_flux(int i,double band){
  if(i < int(observations.size()))
    return observations[i]->get_flux(band);
  else
    return -1;
}

void obs_lib::get_colors(int i,double &c1,double &c2){
  if(i < int(observations.size()))
    observations[i]->get_colors(c1,c2);
  else{
    c1 = -1;
    c2 = -1;
  }
}

void obs_lib::get_all_colors(double* &c1,double* &c2){
  c1 = new double[observations.size()];
  c2 = new double[observations.size()];
  for (unsigned int i=0;i<observations.size();i++){
    observations[i]->get_colors(c1[i],c2[i]);
  }
}

void obs_lib::info(double *b,double *flims,double *ferrs){
  b = new double[3];
  flims = new double[3];
  ferrs = new double[3];

  for(int i=0;i<3;i++){
    b[i] = bands[i];
    flims[i] = flim[i];
    ferrs[i] = ferr[i];
  }
}

obs_lib::~obs_lib(){
  for (unsigned int i=0;i<observations.size();i++){
    delete observations[i];
  }
}
