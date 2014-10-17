#include "sed_lib2.h"

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

sed_lib::sed_lib(string fitsfile){

  color_exp = 0;
  color_zcut=0;
  color_evolution = 1;
  
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
  double temp, scale;
  HDU& header = pInfile->pHDU();

  char buffer[10];
  for (unsigned int i=0;i<lnum;i++) {
    std::string a = "ROW";
    sprintf(buffer,"%i",i+1);
    a.append(buffer);
    try{header.readKey(a,temp);}
    catch(HDU::NoSuchKeyword){
      printf("Keyword %s missing from header (%s)\n",a.c_str(),fitsfile.c_str());
      exit(1);
    }
    lums[i] = temp;
    inds[i] = i;
  }
  
  try{header.readKey("SCALE",temp);}
  catch(HDU::NoSuchKeyword){
    printf("Keyword SCALE missing from header (%s)\n",fitsfile.c_str());
    exit(1);
  }
  scale = pow(10,temp);
  
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline,lnum);
  gsl_spline_init(spline,lums,inds,lnum);
  
  double fluxes[bandnum];
  sed *new_sed;
  
  for (unsigned int fi=0;fi<bandnum;fi++){
    bands[fi]=contents[fi]*scale;
  }
  
  brange[0] = bands[0];
  brange[1] = bands[bandnum-1];
  
  for (unsigned int fj=1;fj<(lnum+1);fj++){
    for (unsigned int fi=0;fi<bandnum;fi++)
      fluxes[fi]=contents[fi+bandnum*fj];
    new_sed = new sed(fluxes,bands,bandnum);
    seds.push_back(new_sed);
  }

  color_exp = 0;

  delete[] bands;
  delete pInfile;
}

//determine the fraction of sources that will be flagged as "AGN" in any given lum/redshift bin
float sed_lib::get_agn_frac(double lum, double redshift, short filter_id)

  float av_lglx,av_lgl6um,av_lgAGNir;
//  float agn_frac;

  //the Chen,Hickox et al. 2013 relation between average(Log(Lx)) and Log(Lir)

  av_lglx=30.37+1.05*lum;
  
  //the Lutz et al. 2004 relation between Log(Lx) and Log(L6um) for Seybert 1s

  av_lgl6um=av_lglx-0.41;

  //the average lgl6 to lg(Lagn_ir) conversion, for now a placeholder only, will need to calculate properly

  av_lgAGNir=av_lgl6um+0.4;

  //the logic here is that if the average Lir_agn is comparable to the Lir of the given bin, essentially 100% of the galaxies are AGN and scales from there.
  return pow(10,(av_lgAGNir-lum));

}
