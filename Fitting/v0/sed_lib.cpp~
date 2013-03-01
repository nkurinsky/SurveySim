#include "sed_lib.h"

sed::sed(){
  fluxes = NULL;
  interp_init = false;
}

sed::sed(double * f,double *bands, int bnum){
  fluxes = new double[bnum];
  for (int i=0;i<bnum;i++){
    fluxes[i] = f[i];

  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline,bnum);
  gsl_spline_init(spline,bands,f,bnum);
  interp_init = true;
}

double sed::get_flux(double band){
  if (interp_init)
    return gsl_spline_eval(spline,band,acc);
  else{
    printf("%s \n","ERROR: SED Model Not Initialized");
    return -1;
  }
}

sed::~sed(){
  if(interp_init){
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
  if (fluxes != NULL)
    delete [] fluxes;
}

//modify to read Anna's type of file
model_lib::sed_lib(string fitsfile){
  
  FITS *pInfile;  
  pInfile = new FITS(fitsfile,Read);

  std::valarray<double> contents;
  PHDU& image = pInfile->pHDU();
  image.read(contents); //store contents of image into the array "seds"

  bandnum = image.axis(0);
  lums.resize(image.axis(1)-1); //since the first row here is the lambda array
  seds.resize(image.axis(1)-1);
  bands = new double[bandnum];

 //read in the total IR luminosities associated with the different z=0 SED templates
  double temp;
  HDU& header = pInfile->pHDU();
  
  for (int i=0;i<nlum;i++) {
    std::string a = "ROW";
    //convert integer to string...may not be best way
    a+= static_cast<ostringstream*>( &(ostringstream() << i+1) )->str();
    header.readKey(a,temp);
    lums[i]=temp;  
  }

  double fluxes[bandnum];
  sed *new_sed;
  
  for (int fi=0;fi<bandnum;fi++){
    bands[fi]=contents[fi];
  }
  
  for (int fj=1;fj<(nlum+1);fj++){
    for (int fi=0;fi<bandnum;fi++)
      fluxes[fi]=contents[fi+bandnum*fj];
    new_sed = new sed(fluxes,bands,bandnum);
    seds[i] = new_sed;
  }
  
  delete pInfile;
  
}

//needs to interpolate to correct luminosity 
//make sure lums sorted in order!!! 
double sed_lib::get_flux(double lum,double band){
  if(lum < seds.size())
    return seds[i]->get_flux(band);
  else
    return -1;
}

sed_lib::~sed_lib(){
  for (unsigned int i=0;i<seds.size();i++){
    delete seds[i];
  }
}
