#include "functions.h"

axis_type set_axis_type(string &keyvalue){
  string saxis(toLower(keyvalue));

  if(saxis == "colorf1f3")
    return ColorF1F3;
  else if(saxis == "colorf2f3")
    return ColorF2F3;
  else if(saxis == "colorf1f2")
    return ColorF1F2;
  else if(saxis == "flux1")
    return Flux1;
  else if(saxis == "flux2")
    return Flux2;
  else if(saxis == "flux3")
    return Flux3;
  else{
    printf("set_axis_type: keyvalue \"%s\" invalid\n",keyvalue.c_str());
    throw(-1);
  }    
}

string get_axis_type(axis_type opt){
  switch(opt){
  case ColorF1F3:
    return "ColorF1F3";
  case ColorF2F3:
    return "ColorF2F3";
  case ColorF1F2:
    return "ColorF1F2";
  case Flux1:
    return "Flux1";
  case Flux2:
    return "Flux2";
  case Flux3:
    return "Flux3";
  default:
    printf("get_axis_type: option %i invalid\n",opt);
  }
  return "";
}

double metric_value(const double& f1,const double &f2,const double &f3,const axis_type &opt){

  switch(opt){
  case ColorF1F3:
    return log10(f3/f1);
  case ColorF2F3:
    return log10(f3/f2);
  case ColorF1F2:
    return log10(f2/f1);
  case Flux1:
    return log10(f1);
  case Flux2:
    return log10(f2);
  case Flux3:
    return log10(f3);
  default:
    printf("Simulator::metric_value Error: Unknown option\n");
    return -99.0;
  }
  
}

string toLower(const string &oldstr){
  string newstr(oldstr);
  for(string::iterator it = newstr.begin(); it != newstr.end(); it++)
    *it = tolower(*it);
  return newstr;
}

Configuration::Configuration(int argc, char *argv[]){
  if(argc < 4){
    printf("ERROR: Invalid number of arguments.\n");
    printf("Calling sequence should be \"%s obsfile modfile sedfile [output]\"\n",argv[0]);
    exit(1);
  }
  
  //File names passed in by Widget
  obsfile = argv[1];
  modfile = argv[2];
  sedfile = argv[3];
  //If outfile specified as argument, change from default
  if(argc > 4)
    outfile = argv[4];
  else
    outfile = "output.fits";
  
  load();
}

void Configuration::print(){

  printf("\n>> Printing Configuration <<\n\n");
  //output found settings
  if(oprint)
    printf("Using Verbose Output\n\n");
  else
    printf("Using Concise Output\n\n");

  printf("File Settings:\n");
  printf("  Main Settings:  %s\n",modfile.c_str());
  printf("  SED library:    %s\n",sedfile.c_str());
  printf("  Observations:   %s\n",obsfile.c_str());

  printf("\nMCMC Settings:\n");
  printf("  Chain Number         : %lu\n", nchain);
  printf("  Starting Temp        : %5.2f\n",tmax);
  printf("  Ideal Accept Pct     : %4.2f\n",idealpct);
  printf("  Burn-in Step         : %lu\n",burn_step);
  printf("  Convergence Step     : %lu\n",conv_step);
  printf("  Run:Burn-In          : %lu\n",burn_ratio);
  printf("  Convergence Criterion: %4.2f\n",rmax);
  printf("  Confidence Interval  : %4.2f\n",a_ci);
 
  printf("\nSimulation Settings:\n");
  printf("  Run Number Max       : %lu\n",runs);
  printf("  Number Redshift Bins : %i\n",nz);
  printf("  Redshift Bin Width   : %f\n",dz);
  printf("  Area [square deg]    : %f\n",area);

  string pnames[] = {"PHI0","L0","ALPHA","BETA","P","Q","ZCUT"};  
  printf("\nLuminosity Function Parameter Settings:\n");
  printf("Parameter\tStart\t Min \t Max \tStep \tFit\n");
  for(int i=0;i<LUMPARS;i++){
    printf("%9s\t%5.2f\t%5.2f\t%5.2f\t%s\n",
	   pnames[i].c_str(),
	   LFParameters[i][value],
	   LFParameters[i][min],
	   LFParameters[i][max],
	   LFParameters[i][fixed] == 0.0 ? "True" : "False"
	   );
  }
  printf("%9s\t%5.2f\t%5.2f\t%5.2f\t%s\n\n",
	 "CEXP",
	 colorEvolution[value],
	 colorEvolution[min],
	 colorEvolution[max],
	 colorEvolution[fixed] == 0.0 ? "True" : "False"
	 );

  printf("\nNumber Unfixed Parameters: %lu\n",nparams);
  
}

double Configuration::areaSteradian() const{
  return area*pow((M_PI/180.0),2.0);
}

void Configuration::LFparameters(double lpars[], short type){
  for(int i=0;i<LUMPARS;i++)
    lpars[i] = LFParameters[i][type];
}

void Configuration::load(){
  
  //IO variables
  string pnames[] = {"PHI0","L0","ALPHA","BETA","P","Q","ZCUT"};  
  string suffix[] = {"","_FIX","_MIN","_MAX"};
  string tag;
  
  axes[0] = ColorF1F3;
  axes[1] = ColorF2F3;
  
  std::unique_ptr<CCfits::FITS> pInfile;
  try{
    pInfile.reset(new CCfits::FITS(modfile,CCfits::Read,true));
  }
  catch(...){
    printf("Error reading FITS file %s\n",modfile.c_str());
    exit(1);
  }

  //reading primary header from fits file
  CCfits::HDU& tab = pInfile->pHDU();

  string stemp;
  try{
    tab.readKey("AXIS1",stemp);
    axes[0] = set_axis_type(stemp);
  }
  catch(CCfits::HDU::NoSuchKeyword){
    printf("Keyword \"Axis1\" not specified, defaulting to ColorF1F3\n");
  }
  catch(...){
    printf("Value of \"Axis1\" invalid, defaulting to ColorF1F3\n");
  }

  try{
    tab.readKey("AXIS2",stemp);
    axes[1] = set_axis_type(stemp);
  }
  catch(CCfits::HDU::NoSuchKeyword){
    printf("Keyword \"Axis2\" not specified, defaulting to ColorF2F3\n");
  }
  catch(...){
    printf("Value of \"Axis2\" invalid, defaulting to ColorF2F3\n");
  }

  double rtemp;
  //Read in simulation settings
  tab.readKey("RUNS",rtemp);
  runs = static_cast<unsigned long>(rtemp);
  tab.readKey("ZMIN",zmin);
  tab.readKey("ZMAX",zmax);
  tab.readKey("DZ",dz);
  tab.readKey("AREA",area);
  //generalize:
  nsim = 1000;

  tab.readKey("NCHAIN",rtemp);
  nchain = static_cast<unsigned long>(rtemp);
  tab.readKey("TMAX",tmax);
  tab.readKey("TSCALE",tscale);
  tab.readKey("ANN_PCT",idealpct);
  tab.readKey("ANN_RNG",annrng);
  tab.readKey("BURN_STE",rtemp);
  burn_step = static_cast<unsigned long>(rtemp);
  tab.readKey("CONV_STE",rtemp);
  conv_step = static_cast<unsigned long>(rtemp);
  tab.readKey("BURNVRUN",rtemp);
  burn_ratio = static_cast<unsigned long>(rtemp);
  tab.readKey("CONV_RMA",rmax);
  tab.readKey("CONV_CON",a_ci);
  tab.readKey("PRINT",rtemp);
  oprint = rtemp == 0.0 ? true : false;

  //compute redshift bin number from FITS values
  nz = (zmax-zmin)/dz;

  //=================================================================  
  //Read-in Luminosity Function Parameters
  //-----------------------------------------------------------------

  for(int i=0;i<LUMPARS;i++){
    for(int j=0;j<4;j++){
      tag = pnames[i]+suffix[j];
      tag = tag.substr(0,8);
      tab.readKey(tag,LFParameters[i][j]);
    }
    if(LFParameters[i][fixed] == 0)
      param_inds.push_back(i);
  }
  
  //color evolution parameters
  for(int j=0;j<4;j++){
    tag = "CEXP"+suffix[j];
    tag = tag.substr(0,8);
    tab.readKey(tag,colorEvolution[j]);
  }

  nparams = param_inds.size();
  cind = nparams;
  if(colorEvolution[fixed] == 0){
    nparams++;
    vary_cexp = true;
  }
  else 
    vary_cexp = false;
  
  burn_num = runs/burn_ratio;
}

RandomNumberGenerator::RandomNumberGenerator(){
  gsl_rng_default_seed=time(NULL);
  T=gsl_rng_default;
  r=gsl_rng_alloc(T);
}

RandomNumberGenerator::~RandomNumberGenerator(){
  gsl_rng_free (r);
}

double RandomNumberGenerator::gaussian(double mean, double sigma, double min, double max){
  static double temp;
  do{
    temp = mean+gsl_ran_gaussian(r,sigma);
  }while((temp < min) or (temp > max));
  
  return temp;
}

double RandomNumberGenerator::flat(double min, double max){
  return gsl_ran_flat(r,min,max);
}

void RandomNumberGenerator::gaussian_mv(const vector<double> &mean, const vector<vector<double> > &covar, const vector<double> &min, const vector<double> &max, vector<double> &result){
  
  /* multivariate normal distribution random number generator */
  /*
   *	n	dimension of the random vetor
   *	mean	vector of means of size n
   *	var	variance matrix of dimension n x n
   *	result	output variable with a sigle random vector normal distribution generation
   */
  int k;
  int n=mean.size();
  gsl_matrix *_covar = gsl_matrix_alloc(covar.size(),covar[0].size());
  gsl_vector *_result = gsl_vector_calloc(mean.size());
  gsl_vector *_mean = gsl_vector_calloc(mean.size());
  result.resize(mean.size());

  for(k=0;k<n;k++){
    for(int j=0;j<n;j++){
      gsl_matrix_set(_covar,k,j,covar[k][j]);
    }
    gsl_vector_set(_mean, k, mean[k]);
  }

  gsl_linalg_cholesky_decomp(_covar);
  
  bool in_range;
  do{
    for(k=0; k<n; k++)
      gsl_vector_set( _result, k, gsl_ran_ugaussian(r) );
    
    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, _covar, _result);
    gsl_vector_add(_result,_mean);
    
    in_range = true;
    for(k=0; k<n; k++){
      if(gsl_vector_get(_result, k) < min[k] or gsl_vector_get(_result, k) > max[k]){
	in_range = false;
	k=n+1;
      }
    }
  }while(not in_range);

    for(k=0; k<n; k++){
      result[k] = gsl_vector_get(_result, k);
    }

  gsl_matrix_free(_covar);
  gsl_vector_free(_result);
  gsl_vector_free(_mean);
  
  return;
}
