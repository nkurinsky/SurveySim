#include "sed_lib.h"

sed::sed(){
  fluxes = NULL;
  interp_init = false;
}

sed::sed(double * f,double *bands){
  for (int i=0;i<3;i++){
    fluxes[i] = f[i];

  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline,3);
  gsl_spline_init(spline,bands,f,3);
  interp_init = true;
}

double sed::get_flux(double band){
  return gsl_spline_eval(spline,band,acc);
}

sed::~sed(){
  if(interp_init){
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
  delete [] fluxes;
}

//modify to read Anna's type of file
model_lib::sed_lib(string fitsfile){
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
  seds.resize(image.axis(1));

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
    
    sed *new_sed;
    double fluxes[3];
    
    for (unsigned int i=0;i<seds.size();i++){
      for (int j=0;j<3;j++)
	fluxes[j]=contents[i*3+j];
      new_sed = new sed(fluxes,bands);
      seds[i] = new_mods;
    }
  }
  
  delete pInfile;
}
  
double sed_lib::get_flux(double lum,double band){
  if(i < int(seds.size()))
    return seds[i]->get_flux(band);
  else
    return -1;
}

sed_lib::~sed_lib(){
  for (unsigned int i=0;i<seds.size();i++){
    delete seds[i];
  }
}
