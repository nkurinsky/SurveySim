#include "model_lib.h"

mods::mods(){
  c1 = -99;
  c2 = -99;
  interp_init = false;
}

mods::mods(double * f,double *bands){
  bool do_colors = true;
  for (int i=0;i<3;i++){
    fluxes[i] = f[i];
    if((bands[i] <= 0.0000000E0) or (f[i] <= 0.00000000E0))
      do_colors = false;
  }

  if(do_colors){
    c1 = log10(fluxes[2]/fluxes[0]) ; //get_color(fluxes[2],fluxes[0],bands[2],bands[0]);
    c2 = log10(fluxes[2]/fluxes[1]) ; //get_color(fluxes[2],fluxes[1],bands[2],bands[1]);
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

double mods::get_flux(double band){
  return gsl_spline_eval(spline,band,acc);
}

void mods::get_colors(double &c1,double &c2){
  c1 = this->c1;
  c2 = this->c2;
}

mods::~mods(){
  if(interp_init){
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
}

model_lib::model_lib(string fitsfile){
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
  models.resize(image.axis(1));

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
    
    mods *new_mods;
    double fluxes[3];
    
    for (unsigned int i=0;i<models.size();i++){
      for (int j=0;j<3;j++)
	fluxes[j]=contents[i*3+j];
      new_mods = new mods(fluxes,bands);
      models[i] = new_mods;
    }
  }
  
  delete pInfile;
}
  
double model_lib::get_flux(int i,double band){
  if(i < int(models.size()))
    return models[i]->get_flux(band);
  else
    return -1;
}

void model_lib::get_colors(int i,double &c1,double &c2){
  if(i < int(models.size()))
    models[i]->get_colors(c1,c2);
  else{
    c1 = -1;
    c2 = -1;
  }
}

void model_lib::get_all_colors(double* &c1,double* &c2){
  c1 = new double[models.size()];
  c2 = new double[models.size()];
  for (unsigned int i=0;i<models.size();i++){
    models[i]->get_colors(c1[i],c2[i]);
  }
}

void model_lib::info(double *b,double *flims,double *ferrs){
  b = new double[3];
  flims = new double[3];
  ferrs = new double[3];

  for(int i=0;i<3;i++){
    b[i] = bands[i];
    flims[i] = flim[i];
    ferrs[i] = ferr[i];
  }
}

model_lib::~model_lib(){
  for (unsigned int i=0;i<models.size();i++){
    delete models[i];
  }
}
