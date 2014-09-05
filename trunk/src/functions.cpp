#include "functions.h"

//Computes spectral color in log space from flux and band inputs
double get_color(const double &f1,const double &f2){
  return log10(f1/f2);
}

double metric_value(const double& f1,const double &f2,const double &f3,const axis_type &opt){

  switch(opt){
  case ColorF1F3:
    return get_color(f3,f1);
  case ColorF2F3:
    return get_color(f3,f2);
  case ColorF1F2:
    return get_color(f2,f1);
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
  
  char * ffile = getenv("FILTERFILE");
  if(ffile != NULL)
    filterfile = static_cast<string>(ffile);
  else
    filterfile = "/usr/local/surveysim/filters/filterlib.txt";

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
  printf("  Filter library: %s\n",filterfile.c_str());

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

double Configuration::areaSteradian(){
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
    pInfile.reset(new CCfits::FITS(modfile,CCfits::Read));
  }
  catch(...){
    printf("Error reading FITS file %s\n",modfile.c_str());
    exit(1);
  }
  
  //reading primary header from fits file
  CCfits::HDU& tab = pInfile->pHDU();

  double rtemp;
  //Read in simulation settings
  tab.readKey("RUNS",rtemp);
  runs = static_cast<unsigned long>(rtemp);
  tab.readKey("ZMIN",zmin);
  tab.readKey("ZMAX",zmax);
  tab.readKey("DZ",dz);
  tab.readKey("AREA",area);

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
