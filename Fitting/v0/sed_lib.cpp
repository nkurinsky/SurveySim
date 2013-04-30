#include "sed_lib.h"

sed::sed(){
  fluxes = NULL;
  interp_init = false;
}

sed::sed(double * f,double *bands, int bnum){
  fluxes = new double[bnum];
  for (int i=0;i<bnum;i++){
    fluxes[i] = f[i];}

  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline,bnum);
  gsl_spline_init(spline,bands,fluxes,bnum);
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
sed_lib::sed_lib(string fitsfile){
  
  double *bands;
  FITS *pInfile;  
  pInfile = new FITS(fitsfile,Read);

  std::valarray<double> contents;
  PHDU& image = pInfile->pHDU();
  image.read(contents); //store contents of image into the array "seds"

  bandnum = image.axis(0);
  lnum = image.axis(1)-1; //since the first row here is the lambda array
  lums = new double[lnum];
  double inds[lnum];
  seds.reserve(lnum);
  bands = new double[bandnum];

 //read in the total IR luminosities associated with the different z=0 SED templates
  double temp;
  HDU& header = pInfile->pHDU();
  
  printf("%s\n","Initialize Luminosities");
  for (unsigned int i=0;i<lnum;i++) {
    std::string a = "ROW";
    //convert integer to string...may not be best way
    a+= static_cast<ostringstream*>( &(ostringstream() << i+1) )->str();
    header.readKey(a,temp);
    lums[i] = temp;
    inds[i] = i;
  }
  
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline,lnum);
  gsl_spline_init(spline,lums,inds,lnum);

  double fluxes[bandnum];
  sed *new_sed;
  
  for (unsigned int fi=0;fi<bandnum;fi++){
    bands[fi]=contents[fi];
  }
  
  brange[0] = bands[0];
  brange[1] = bands[bandnum-1];
  
  for (unsigned int fj=1;fj<(lnum+1);fj++){
    for (unsigned int fi=0;fi<bandnum;fi++)
      fluxes[fi]=contents[fi+bandnum*fj];
    new_sed = new sed(fluxes,bands,bandnum);
    seds.push_back(new_sed);
  }

  delete[] bands;
  delete pInfile;
}

//needs to interpolate to correct luminosity 
//rounds to nearest template, in future may average/interp
//make sure lums sorted in order!!! currently an assumption
double sed_lib::get_flux(double lum,double band){
  static int i = 0;
  i = int(floor(0.5+ gsl_spline_eval(spline,lum,acc)));
  if (i < 0)
    i = 0;
  if (i >= int(lnum))
    i = int(lnum)-1;
  if (seds[i] != NULL){
    if((band >= brange[0]) and (band <= brange[1]))
      return seds[i]->get_flux(band);
    else{
      printf("%s\n%s\n","ERROR: Band out of model range.","Check that the obs bands are within range and unit consistent");
      return -1;
    }
  }
  else{
    printf("%s %i \n","ERROR: Unitialized Model Called!!!",i);
    return -1;
  }
}

void sed_lib::get_lums(double luminosities[]){
  for (unsigned int i=0;i < lnum; i++)
    luminosities[i] = lums[i];
}

sed_lib::~sed_lib(){
  for (unsigned int i=0;i<seds.size();i++){
    delete seds[i];
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  delete[] lums;
}
