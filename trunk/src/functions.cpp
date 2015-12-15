#include "functions.h"

colsel_type set_colsel_type(string &keyvalue){
  //string scolsel(toLower(keyvalue));
  string scolsel=keyvalue;
  if(scolsel == "None")
    return None;
  else if(scolsel == "mag1_mag2")
    return mag1_mag2;
  else if(scolsel == "ColF1F2")
    return ColF1F2;
  else if(scolsel == "ColF1F3")
    return ColF1F3;
  else if(scolsel == "ColF2F3")
    return ColF2F3;
  else if(scolsel == "ColF2F1")
    return ColF2F1;
  else{
    printf("set_colsel_type: keyvalue \"%s\" invalid\n",keyvalue.c_str());
  throw(-1);
  }
}

axis_type set_axis_type(string &keyvalue){
  string saxis(toLower(keyvalue));

  if(saxis == "colorf1f3")
    return ColorF1F3;
  else if(saxis == "colorf2f3")
    return ColorF2F3;
  else if(saxis == "colorf1f2")
    return ColorF1F2;
  else if(saxis == "colorf2f1")
    return ColorF2F1;
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
  case ColorF2F1:
    return "ColorF2F1";
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

string get_colsel_type(colsel_type  opt){
  switch(opt){
  case None:
    return "None";
  case mag1_mag2:
    return "mag1_mag2";
  case ColF1F2:
    return "ColF1F2";
  case ColF1F3:
    return "ColF1F3";
  case ColF2F3:
    return "ColF2F3";
  case ColF2F1:
    return "ColF2F1";
  default:
    printf("get_colsel_type: option %i invalid\n",opt);
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
  case ColorF2F1:
    return log10(f1/f2);
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

  if((argc > 1) and (argv[1][0] == '-') and (argv[1][1] == 'h')){
    printf("Calling sequence:\"%s modfile sedfile [obsfile] [-o outfile] [-v]\"\n",argv[0]);
    printf("\tmodfile - Model file containing parameters and filters\n");
    printf("\tsedfile - SED template file\n");
    printf("\tobsfile - Observation file, table of fluxes\n");
    printf("\t-o outfile - Name of output file, defaults to \"output.fits\"\n");
    printf("\t-v - Force verbose output\n");
    exit(0);
  }
  
  if(argc < 3){
    printf("ERROR: Invalid number of arguments.\n");
    printf("Calling sequence should be \"%s modfile sedfile [obsfile] [-o output] [-v]\"\n",argv[0]);
    exit(1);
  }
  
  //File names passed in by Widget
  modfile = argv[1];
  sedfile = argv[2];
  outfile="output.fits";
  //If outfile specified as argument, change from default
  int i=3;

  simflag=true;
  bool force_verbose=false;
  while(i<argc){
    if((argv[i][0] == '-') and (argv[i][1] == 'v')){
      force_verbose=true;
      i++;
    }
    else if((argv[i][0] == '-') and (argv[i][1] == 'o')){
      if(i+1 < argc){
	outfile = argv[i+1];
	i+=2;
      }
      else{
	fprintf(stderr,"ERROR: Please enter output file name after -o\n");
	exit(1);
      }
    }
    else if (i >= 3){
      simflag=false;
      obsfile=argv[i];
      i++;
    }
    else{
      fprintf(stderr,"ERROR: Unkown argument \"%s\" after obsfile argument\n",argv[i]);
      exit(1);
    }
  }
  
  load();

  if(force_verbose)
    oprint=2;
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
  printf("  Starting Temp        : %5.2f\n",temp);
  printf("  Learning Rate        : %5.2f\n",learningRate);
  printf("  Ideal Accept Pct     : %4.2f\n",idealpct);
  printf("  Burn-in Step         : %lu\n",burn_step);
  printf("  Convergence Step     : %lu\n",conv_step);
  printf("  Run:Burn-In          : %lu\n",burn_ratio);
  printf("  Convergence Criterion: %4.2f\n",rmax);
  printf("  Confidence Interval  : %4.2f%%\n",1-a_ci);
 
  printf("\nSimulation Settings:\n");
  printf("  Run Number Max       : %lu\n",runs);
  printf("  Number Redshift Bins : %i\n",nz);
  printf("  Redshift Bin Width   : %f\n",dz);
  printf("  Area [square deg]    : %f\n",area);

  printf("\nLuminosity Function Form: ");
  switch(lfDist){
  case LF::distribution::Schecter:
    printf("Schecter Function");
    break;
  case LF::distribution::DoublePowerLaw:
    printf("Double Power Law");
    break;
  case LF::distribution::ModifiedSchecter:
    printf("Modified Schecter Function");
    break;
  }
  printf("\n");

  string pnames[] = {"PHI0","L0","ALPHA","BETA","P","Q","P2","Q2","ZBP","ZBQ","FA0","T1","T2","ZBT","FCOMP","FCOLD"};  
  printf("\nLuminosity Function Parameter Settings:\n");
  printf("Parameter\tStart\t Min \t Max \tFit\n");
  for(int i=0;i<LUMPARS;i++){
    printf("%9s\t%5.2f\t%5.2f\t%5.2f\t%s\n",
	   pnames[i].c_str(),
	   LFParameters[i][value],
	   LFParameters[i][min],
	   LFParameters[i][max],
	   LFParameters[i][fixed] == 0.0 ? "True" : "False"
	   );
  }
  printf("%9s\t%5.2f\t%5.2f\t%5.2f\t%s\n",
  	 "CEXP",
  	 colorEvolution[value],
  	 colorEvolution[min],
  	 colorEvolution[max],
  	 colorEvolution[fixed] == 0.0 ? "True" : "False"
	 );
  printf("%9s\t%5.2f\t%5.2f\t%5.2f\t%s\n\n",
  	 "ZBC",
  	 colorZCut[value],
  	 colorZCut[min],
  	 colorZCut[max],
  	 colorZCut[fixed] == 0.0 ? "True" : "False"
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
  string pnames[] = {"PHI0","L0","ALPHA","BETA","P","Q","P2","Q2","ZBP","ZBQ","FA0","T1","T2","ZBT","FCOMP","FCOLD"};  
  string suffix[] = {"","_FIX","_MIN","_MAX"};
  string tag;
  
  axes[0] = ColorF1F3;
  axes[1] = ColorF2F3;
  lfDist = LF::DoublePowerLaw;

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

  try{
    tab.readKey("COLSEL",stemp);
    colsel = set_colsel_type(stemp);
  }
  catch(CCfits::HDU::NoSuchKeyword){
    printf("Keyword \"COLSEL\" not specified, defaulting to None\n");
  }
  catch(...){
    printf("Value of \"COLSEL\" invalid, defaulting to None\n");
  }

  double rtemp;
  //Read in simulation settings
  try{
    tab.readKey("RUNS",rtemp);
    runs = static_cast<unsigned long>(rtemp);
  }
  catch(CCfits::HDU::NoSuchKeyword){
    printf("Value of \"RUNS\" invalud, defaulting to 1000\n");
    runs=1000;
  }
  try{
    tab.readKey("SIGSIZE",sigmaSize);
  }
  catch(CCfits::HDU::NoSuchKeyword){
    sigmaSize=12.0;
  }

  try{
    tab.readKey("ZMIN",zmin);
    tab.readKey("ZMAX",zmax);
    tab.readKey("DZ",dz);
    tab.readKey("AREA",area);
    //generalize:
    nsim = 1000;
    
    tab.readKey("NCHAIN",rtemp);
    nchain = static_cast<unsigned long>(rtemp);
    tab.readKey("TEMP",temp);
    tab.readKey("LRATE",learningRate);
    tab.readKey("ANN_PCT",idealpct);
    tab.readKey("ANN_RNG",annrng);
    tab.readKey("BURN_STE",rtemp);
    burn_step = static_cast<unsigned long>(rtemp);
    tab.readKey("CONV_STE",rtemp);
    conv_step = static_cast<unsigned long>(rtemp);
    tab.readKey("BURNVRUN",rtemp);
    burn_ratio = static_cast<unsigned long>(rtemp);
    tab.readKey("CONV_CON",rtemp);
    a_ci=1-rtemp;
    rmax=1+a_ci;
    tab.readKey("PRINT",rtemp);
    oprint = rtemp; // == 0.0 ? true : false; 

    tab.readKey("LF_FORM",rtemp);
    switch(static_cast<int>(rtemp)){
    case 0:
      lfDist=LF::DoublePowerLaw;
      break;
    case 1:
      lfDist=LF::ModifiedSchecter;
      break;
    case 2:
      lfDist=LF::Schecter;
      break;
    default:
      cerr << "Error: Invalid LF Form " << rtemp << endl;
      exit(90);
    }

  }
  catch(CCfits::FitsException e){
    printf("Error reading model file\nPlease check that all keywords are present\n");
    cerr << e.message() << endl;
    exit(10);
  }

  //compute redshift bin number from FITS values
  nz = (zmax-zmin)/dz;

  //=================================================================  
  //Read-in Luminosity Function Parameters
  //-----------------------------------------------------------------

  for(int i=0;i<LUMPARS;i++){
    for(int j=0;j<4;j++){
      tag = pnames[i]+suffix[j];
      tag = tag.substr(0,8);
      try{
	tab.readKey(tag,LFParameters[i][j]);
      }
      catch(CCfits::HDU::NoSuchKeyword){
	printf("Error reading keyword %s\n",tag.c_str());
	exit(12);
      }
    }
    if(LFParameters[i][fixed] == 0)
      param_inds.push_back(i);
  }
  
  //color evolution parameters
  for(int j=0;j<4;j++){
    tag = "CEXP"+suffix[j];
    tag = tag.substr(0,8);
    try{
      tab.readKey(tag,colorEvolution[j]);
    }
    catch(CCfits::HDU::NoSuchKeyword){
      printf("Error reading keyword %s\n",tag.c_str());
      exit(12);
    }
  }

  for(int j=0;j<4;j++){
    tag = "ZBC"+suffix[j];
    tag = tag.substr(0,8);
    try{
      tab.readKey(tag,colorZCut[j]);
    }
    catch(CCfits::HDU::NoSuchKeyword){
      printf("Error reading keyword %s\n",tag.c_str());
      exit(12);
    }
  }

  nparams = param_inds.size();
  cind = nparams;
  if(colorEvolution[fixed] == 0){
    nparams++;
    vary_cexp = true;
  }
  else 
    vary_cexp = false;

  zbcind = nparams;
  if(colorZCut[fixed] == 0){
    nparams++;
    vary_zbc = true;
  }
  else 
    vary_zbc = false;
  
  burn_num = runs/burn_ratio;


  //completeness information
  char cnum[2];
  string num;
  for (int i=0;i<3;i++){
    sprintf(cnum,"%i",i+1);
    num = cnum;
    try{tab.readKey("COMP"+num+"N",completeness_n[i]);}
    catch(CCfits::HDU::NoSuchKeyword){
      completeness_n[i]=-1.0;
    }
    try{tab.readKey("COMP"+num+"B",completeness_b[i]);}
    catch(CCfits::HDU::NoSuchKeyword){
      completeness_b[i]=-1.0;
    }
    try{tab.readKey("COMP"+num+"M",completeness_m[i]);}
    catch(CCfits::HDU::NoSuchKeyword){
      completeness_m[i]=0.0;
    }
  }

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

  temp = mean+gsl_ran_gaussian(r,sigma);
  if((temp >= min) and (temp <= max))
    return temp;
  else if(temp > max)
    return max;
  else
    return min;
}

double RandomNumberGenerator::poisson(double lambda){
  static double temp;

  temp = lambda+gsl_ran_gaussian(r,sqrt(lambda));
  if(temp >= 0)
    return temp;
  else
    return 0;
}

double RandomNumberGenerator::flat(double min, double max){
  return gsl_ran_flat(r,min,max);
}

double RandomNumberGenerator::triangular(double min, double max, double fmin, double fmax){
  double u = gsl_ran_flat(r,0,1);
  double r = (fmax/fmin);
  double slope = (sqrt(1.0+(pow(r,2.0)-1.0)*u)-1.0)/(r-1.0);

  return slope*(max-min)+min;
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

CompletenessCurve::CompletenessCurve(double n, double m, double b){
  static int bnum=1;
  //invalid to have -n
  if(n > 0.0){
    _n=n;
    complete=false;
  }
  else{
    _n=0.0;
    complete=true;
  }
  
  //m has infinite domain
  _M=m;

  //doesn't make sense to have -b
  if(b > 0.0)
    _B=b;
  else
    _B=1e-10;
  
  scaled=false;
  if(not complete){
    _ulim=_M-_B*log(pow(0.999,-1.0/_n)-1);
    _fScale=50.0/_ulim;
    if(_fScale < 1.0)
      _fScale=1.0;   
    else
      scaled=true;
  }
  else{
    _ulim=0;
    _fScale=1.0;
  }

  cout << "Completeness (" << bnum << "): " << _B << " " << _M << " " << _n << ", Max= " << _ulim << ", Scale=" << 1.0/_fScale << endl;
  if(scaled){
    _M*=_fScale;
    _B*=_fScale;
  }
  bnum++;
}

bool CompletenessCurve::accept(double flux){
  
  if(complete){
    return true;
  }

  if(flux > _ulim)
    return true;

  double tFlux;
  if(scaled)
    tFlux=round(flux*_fScale);
  else
    tFlux=round(flux);  

  if(valueStore.count(tFlux) == 0){
    valueStore[tFlux]=pow(1.0+exp((_M-tFlux)/_B),-_n);
  }
  return valueStore[tFlux] >= crng.flat(0.0,1.0);
}
