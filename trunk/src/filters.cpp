/*
  Comment
*/

#include "filters.h"

filter::filter(){
  name="NULL";
  init=false;
  filter_size=0;
  lambda=NULL;
  response=NULL;
  acc=NULL;
  spline=NULL;
}

filter::filter(string filtername, vector<double> band, vector<double> transmission,int logflag){
  name="NULL";
  init=false;
  filter_size=0;
  lambda=NULL;
  response=NULL;
  acc=NULL;
  spline=NULL;
  load(filtername,band,transmission,logflag);
}

bool filter::load(string filtername, vector<double> band, vector<double> transmission, int logflag){
  
  double width = 0;
  
  if(band.size() < 2){
    LOG_CRITICAL(printf("Error: Not enough valid lines in filter file, filter not initialized\n"));
    return false;
  }  
  if(band.size() != transmission.size()){
    LOG_CRITICAL(printf("Error: Filter file improperly formatted (different number of wavelength and response values), filter not initialized\n"));
    return false;
  }
  
  init=true;
  name = filtername;
  filter_size = band.size();
  filter_limits[0] = band.front();
  filter_limits[1] = band.back();

  if(lambda != NULL)
    delete [] lambda;
  if(response != NULL)
    delete [] response;
  
  //convert vectors to arrays for gsl functions and storage
  lambda = new double[filter_size];
  response = new double[filter_size];
  //compute normalization
  width = trap_integrate(band,transmission);
  for(unsigned long i=0;i<filter_size;i++){
    lambda[i] = band[i];
    response[i] = transmission[i]/width;
  }
    
  if(init){
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }   

  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline,filter_size);
  gsl_spline_init(spline,lambda,response,filter_size);
  
  return true;
}

double filter::transmission(double wavelength){
  if (init)
    if((wavelength < low()) or (wavelength > high()))
      return 0;
    else
      return gsl_spline_eval(spline,wavelength,acc);
  else{
    printf("ERROR: Filter Not Initialized\n");
    return -1;
  }
}

double filter::low(){
  if(init)
    return filter_limits[0];
  else
    return 0;
}

double filter::high(){
  if(init)
    return filter_limits[1];
  else
    return -1;
}

void filter::print(bool all){
  if(init){
    printf("%s: Range = %g - %g, Points = %lu\n",name.c_str(),filter_limits[0],filter_limits[1],filter_size);
    if(all){
      printf("Printing all filter data points\n\tWavelength\tResponse\n");
      for(unsigned long i=0;i<filter_size;i++)
	printf("\t%g\t%f\n",lambda[i],response[i]);
    }
  }
  else
    printf("ERROR: Filter not initialized\n");
}

filter::~filter(){
  if(init){
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
  if (lambda != NULL){
    delete [] lambda;
    delete [] response;
  }
}

double filter::trap_integrate(vector<double> lambda, vector<double> response){
  double tsum = 0;
  for(unsigned int i=1;i<lambda.size();i++){
    tsum += (lambda[i] - lambda[i-1])*(response[i]+response[i-1]);
  }
  return tsum/2.0;
}

filter_lib::filter_lib(){
  initialized = false;
}

filter_lib::filter_lib(string fitsfile,int logflag){
  if(load_filters(fitsfile,logflag)){
    initialized = true;
  }
  else{
    printf("ERROR: file \"%s\" not read successfully, library not initialized\n",fitsfile.c_str());
    initialized = false;
  }
}

bool filter_lib::load_filters(string fitsfile,int logflag){

  std::unique_ptr<CCfits::FITS> pInfile;
  try{
    pInfile.reset(new CCfits::FITS(fitsfile,CCfits::Read,true));
  }
  catch(...){
    printf("Error reading FITS file %s\n",fitsfile.c_str());
    exit(1);
  }

  CCfits::HDU& head = pInfile->pHDU();
  CCfits::ExtHDU& filters = pInfile->extension(1);
  int ntcols = filters.numCols();
  if(ntcols != 6){
    printf("Wrong number of filters included (need 3):");
    exit(1);
  }

  vector<double> band;
  vector<double> transmission;
  double units[] = {1.0,1.0,1.0};
  long tablelength = filters.rows();
  
  string bands[] = {"BAND_1","BAND_2","BAND_3"};
  string unitKeys[] = {"UNITS1","UNITS2","UNITS3"};
  string limitKeys[] = {"LIMIT1","LIMIT2","LIMIT3"};
  string errorKeys[] = {"ERROR1","ERROR2","ERROR3"};
  string serrorKeys[] = {"SERROR1","SERROR2","SERROR3"};
  for(int i=0;i<3;i++){

    //reading units, defaults to units in mJy
    try{
      string stemp;
      head.readKey(unitKeys[i],stemp);
      if(toLower(stemp) == "jy"){
	units[i]=1e3;
      }
      else if(toLower(stemp) == "mjy"){
	units[i]=1.0;
      }
      else if(toLower(stemp) == "ujy"){
	units[i]=1e-3;
      }
      else{
	LOG_CRITICAL(printf("Keyword \"%s\" contains unknown unit, defaulting to standard units mJy\n",unitKeys[i].c_str()));
      }
    }
    catch(...){
      LOG_CRITICAL(printf("Error reading keyword \"%s\", defaulting to standard units mJy\n",unitKeys[i].c_str()));
    }

    //reading flux limits, converting units
    try{
      head.readKey(limitKeys[i],limits[i]);
      if(limits[i]<0){
      	LOG_CRITICAL(printf("Warning: keyword \"%s\" contains negative value, defaulting to standard value 0\n",limitKeys[i].c_str()));
      limits[i]=0.0;
      }
      else
	limits[i]*=units[i];
    }
    catch(...){
      printf("Error reading keyword \"%s\"\n",limitKeys[i].c_str());
      exit(10);
    }
    
    //reading errors, converting units
    try{
      head.readKey(errorKeys[i],errors[i]);
      //if(errors[i]<0){
      //	LOG_CRITICAL(printf("Warning: keyword \"%s\" contains negative value, defaulting to standard value 0\n",errorKeys[i].c_str()));
      // errors[i]=0.0;
      //}
      //else
	errors[i]*=units[i];
    }
    catch(...){
      printf("Error reading keyword \"%s\"\n",errorKeys[i].c_str());
      exit(10);
    }

    //reading skew errors, converting units
    try{
      head.readKey(serrorKeys[i],skew_errors[i]);
      if(skew_errors[i]<0){
	LOG_CRITICAL(printf("Warning: keyword \"%s\" contains negative value, defaulting to standard value 0\n",serrorKeys[i].c_str()));
	skew_errors[i]=0.0;
      }
      else
	skew_errors[i]*=units[i];
    }
    catch(...){
      printf("Error reading keyword \"%s\"\n",serrorKeys[i].c_str());
      exit(10);
    }
  }

  double scale = -10;
  try{
    double dtemp;
    filters.readKey("LSCALE", dtemp);
    scale = dtemp;
  }
  catch(...){
    LOG_CRITICAL(printf("Error reading keyword \"LSCALE\", defaulting to Angstroms\n"));
  }
  double exp_scale = pow(10,scale);

  double nullvalue=-1;
  for(int i=1;i<7;i+=2){
    filters.column(i).read(band,1,tablelength,&nullvalue);
    filters.column(i+1).read(transmission,1,tablelength);
    int num = floor(i/2);
    if(band.back() <= 0.0){
      while(band.back() <= 0.0){
	band.pop_back();
	transmission.pop_back();
      }
    }
    for(unsigned int j=0;j<band.size();j++){
      band[j]*=exp_scale;
    }
    LOG_INFO(printf("Loaded Filter %i (Lambda: %8.2e m -> %8.2e m, limit: %5.2e, error: %5.2e, skew error: %5.2e)\n",num+1,band.front(),band.back(),limits[num],errors[num],skew_errors[num]));
    if(not this->filters[num].load(bands[num],band,transmission,logflag) ){
      printf("Error loading filter %i from %s, exiting\n",i,fitsfile.c_str());
      exit(1);
    }
  }

  initialized = true;
  return true;
}

void filter_lib::filter_info(string names[], double fluxLimits[], double fluxErrors[], double skewErrors[]){
  for(int i=0;i<3;i++){
    names[i]=filters[i].get_name();
    fluxLimits[i]=limits[i];
    fluxErrors[i]=errors[i];
    skewErrors[i]=skew_errors[i];
  }
}

filter& filter_lib::get(short num){
  if ((num >=0) and (num <3)){
    return filters[num];
  }
  else{
    LOG_CRITICAL(printf("ERROR: Invalid filter number (%i)\n",num));
    return dummy;
  }
}
