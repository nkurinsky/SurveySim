#include "sed_lib.h"

sed::sed(){
  interp_init = false;
  zdim=false;
}

sed::sed(double * f,double *bands, int bnum) : sed(){
  
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline,bnum);
  gsl_spline_init(spline,bands,f,bnum);
  
  brange[0] = f[0];
  brange[1] = f[bnum-1];
  
  interp_init = true;
}

sed::sed(double * f,double *bands, int bnum, double *z, int znum) : sed(){
  
  alglib::real_1d_array bs;
  alglib::real_1d_array zs;
  alglib::real_1d_array fs;
  
  bs.setcontent(bnum,bands);
  zs.setcontent(znum,z);
  fs.setcontent(bnum*znum,f);

  try{
  alglib::spline2dbuildbilinearv(bs,bnum,zs,znum,fs,1,s);
  }catch(alglib::ap_error e){
    printf("ERROR: Failed to create interpolation for SED\n");
    printf("error msg: %s\n", e.msg.c_str());
    exit(1);
  }  

  brange[0] = bands[0];
  brange[1] = bands[bnum-1];
  zrange[0] = z[0];
  zrange[1] = z[znum-1];
  
  zdim=true;
  interp_init = true;
}

double sed::get_flux(double band, double redshift){
  if (interp_init){
    double x_em = (band/(1+redshift));
    double z = redshift;
    if ((x_em < brange[0]) or (x_em > brange[1])){
      printf("ERROR: emitted lambda %f out of range %f -> %f\n",x_em,brange[0],brange[1]);
      exit(1);
    }
    if (z < zrange[0]){
      z = zrange[0];
    }
    else if(z > zrange[1]){
      z = zrange[1];
    }
    
    if(zdim){
      return alglib::spline2dcalc(s,x_em,z);
    }
    else{
      return gsl_spline_eval(spline,x_em,acc);
    }
  }
  else{
    printf("%s \n","ERROR: SED Model Not Initialized");
    return -1;
  }
}

sed::~sed(){
  if(interp_init and not zdim){
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
}

sed_lib::sed_lib(string fitsfile){

  redshiftq=false;
  typeq=false;

  color_exp = 0;
  color_zcut = 0;
  color_evolution = 1;
  
  double *bands;
  FITS *pInfile;  
  pInfile = new FITS(fitsfile,Read);

  std::valarray<double> contents;
  PHDU& image = pInfile->pHDU();
  image.read(contents); //store contents of image into the array "seds"

  bandnum = image.axis(0);
  lnum = image.axis(1)-1; //since the first row here is the lambda array
  int axnum = image.axes();

  znum=tnum=1;
  if(axnum >= 3){
    redshiftq=true;
    znum = image.axis(2);
    printf("Detected SED redshift dimension, length %i\n",znum);
    if(axnum >= 4){
      typeq=true;
      tnum = image.axis(3);
      printf("Detected SED type dimension, length %i\n",tnum);
    }
  }
  else{
    printf("Redshift independent templates, will use default color evolution\n");
  }

  lums = new double[lnum];
  double inds[lnum];
  seds.reserve(lnum*tnum);
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
  
  double fluxes[bandnum*znum];
  sed *new_sed;
  
  for (unsigned int fi=0;fi<bandnum;fi++){
    bands[fi]=contents[fi]*scale;
  }
  
  brange[0] = bands[0];
  brange[1] = bands[bandnum-1];
  
  double zs[znum];
  if(znum > 1){
    for (unsigned int i=0;i<znum;i++) {
      std::string a = "Z";
      sprintf(buffer,"%i",i);
      a.append(buffer);
      try{header.readKey(a,temp);}
      catch(HDU::NoSuchKeyword){
	printf("Keyword %s missing from header (%s)\n",a.c_str(),fitsfile.c_str());
	exit(1);
      }
      zs[i] = temp;
    }
  }
  else
    zs[0]=0.0;

  for(unsigned int ft=0; ft < tnum; ft++){
    for (unsigned int fj=1;fj<(lnum+1);fj++){
      for(unsigned int fz=0; fz < znum; fz++){
	for (unsigned int fi=0;fi<bandnum;fi++){
	  fluxes[fi+fz*bandnum]=contents[fi+bandnum*(fj+znum*(fz+tnum*ft))];
	}
      }
      if(znum > 1)
	new_sed = new sed(fluxes,bands,bandnum,zs,znum);
      else
	new_sed = new sed(fluxes,bands,bandnum);
      seds.push_back(new_sed);
    }
  }

  color_exp = 0;

  delete[] bands;
  delete pInfile;
}

bool sed_lib::load_filters(string file){
  
  if(filters.load_filters(file)){
    w = gsl_integration_workspace_alloc(SL_INT_SIZE);
    return true;
  }
  
  return false;
}

//needs to interpolate to correct luminosity 
//rounds to nearest template, in future may average/interp
//make sure lums sorted in order!!! currently an assumption
double sed_lib::get_flux(double lum, double band, short sedtype, double redshift){
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
    }
  }
  else{
    printf("%s %i \n","ERROR: Unitialized Model Called!!!",i);
  }
  
  return -1;
}

double sed_lib::get_filter_flux(double lum, double redshift, short sedtype, short filter_id){

  static int i = 0;

  if((sedtype < 0) or (sedtype >= tnum)){
    printf("ERROR: Invalid sedtype \"%i\" out of possible range (0,%i)\n",sedtype,tnum);
    return -1;
  }

  //determine which model to use; pick that closer to 'lum' called
  i = int(floor(0.5+ gsl_spline_eval(spline,lum,acc)));
  
  if (i < 0)
    i = 0;
  else if (i >= int(lnum))
    i = int(lnum)-1;
  
  if (seds[i] != NULL){
    if (filters.init()){
      if((filter_id >= 0) and (filter_id < 3) and (filters.get(filter_id).low() < filters.get(filter_id).high())){
	if((filters.get(filter_id).low() >= brange[0]) and (filters.get(filter_id).high() <= brange[1]))
	  return convolve_filter(i,redshift,sedtype,filter_id);
	else{
	  printf("%s\n%s\n","ERROR: Filter out of model range.","Check that the obs bands are within range and unit consistent");
	}
      }
      else{
	printf("%s %i %s\n","ERROR: Filter number",filter_id,"invalid, 0-2 acceptable, or filter unititialized.");
      }
    }
    else{
      printf("%s %i \n","ERROR: Filter library not initialized!!!",i);
    }
  }
  else{
    printf("%s %i \n","ERROR: Unitialized Model Called!!!",i);
  } 

  return -1;
}

void sed_lib::set_color_evolution(double exp, double zcut){
  if(zcut >= 0)
    color_zcut = zcut;
  color_exp = exp;
  color_evolution = pow(1+zcut,color_exp);
}

double sed_lib::get_dl(){
  static double res;
  if ((lums != NULL) and (lnum >= 2)){
    res = lums[1]-lums[0];
    return (res >= 0) ? res : -1*res;
  }
  else{
    printf("%s\n%s\n","Error: DL requested but luminosities not intialized","Check that sed library valid");
    return 1;
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
  
  if(filters.init())
    gsl_integration_workspace_free(w);
  
  delete[] lums;
}

double sed_lib::convolve_filter(short lum_id, double redshift, short sedtype, short filter_id){

  static map<tuple<short,short,double>,double > fluxes[3];
  tuple<short,short,double> params (lum_id,sedtype,redshift);
  
  //we use a map here to ensure that a given redshift/SED/filter combination is only convolved once; we recompute color evolution every time as this may change in fitting, while filters/SEDs are static.
  if(fluxes[filter_id].count(params) == 0){
    double scale,result,error;
    double bounds[2];
    gsl_function F;
    flux_yield_params p(seds[lum_id+lnum*sedtype],&filters.get(filter_id),redshift);
    F.function = &flux_yield;
    F.params = &p;

    scale = (1.0+redshift)*Wm2Hz_TO_mJy/(4.0*M_PI*pow(lumdist(redshift)*MPC_TO_METER,2.0));

    result = error = -1;
    bounds[0] = filters.get(filter_id).low();
    bounds[1] = filters.get(filter_id).high();
    
    gsl_integration_qags (&F, bounds[0], bounds[1], 0, SL_INT_PRECISION, SL_INT_SIZE, w, &result, &error); 

    result*=scale;
    fluxes[filter_id][params] = result;
    printf("-\n");
  }
  
  if(redshift < color_zcut)
    return fluxes[filter_id][params]*pow((1.0+redshift),color_exp);
  
  return fluxes[filter_id][params]*color_evolution;
}

double flux_yield (double x, void * params) {
  flux_yield_params *p = (flux_yield_params*) params;
  double z = (p->z);
  double flux = (p->SED->get_flux(x,z));
  double transmission = (p->FILT->transmission(x));
  
  return flux*transmission;
}

flux_yield_params::flux_yield_params(sed *mySED, filter *myFILT, double redshift){
  if(mySED == NULL or myFILT == NULL){
    printf("ERROR: Null pointer passed to flux_yield_params constructor\n");
    exit(1);
  }
  
  SED = mySED;
  FILT = myFILT;
  z = redshift;
}
