#include "simulator.h"

simulator::simulator(){
  distributions = NULL;
  luminosity_function = NULL;
  models = NULL;
  observations = NULL;
  diagnostic = NULL;
  distribution_size = 1000;

  redshift_range[0] = 0.0;
  redshift_range[1] = 10.0;

  for (int i=0;i<3;i++){
    bands[i] = -1;
    band_errs[i] = -1;
    flux_limits[i]= -1;
  }
}

simulator::simulator(double b[],double b_err[],double f_lims[]){
  distributions = NULL;
  diagnostic = NULL;
  luminosity_function = NULL;
  models = NULL;
  observations = NULL;
  distribution_size = 1000;

  redshift_range[0] = 0.0;
  redshift_range[1] = 10.0;

  for (int i=0;i<3;i++){
    bands[i] = b[i];
    band_errs[i] = b_err[i];
    flux_limits[i]= f_lims[i];
  }
}

simulator::simulator(lumfunct *lf,string modfile,string obsfile){
  distributions = NULL;
  diagnostic = NULL;
  distribution_size = 1000;

  redshift_range[0] = 0.0;
  redshift_range[1] = 10.0;

  if(lf != NULL)
    luminosity_function = lf;
  else
    cout << "ERROR: NULL Pointer Passed to Simulator" << endl;
  
  models = new model_lib(modfile);
  observations = new obs_lib(obsfile);

  for (int i=0;i<3;i++){
    bands[i] = -1;
    band_errs[i] = -1;
    flux_limits[i]= -1;
  }
}

simulator::simulator(double b[],double b_err[],double f_lims[],lumfunct *lf,string modfile,string obsfile){
  distributions = NULL;
  diagnostic = NULL;
  distribution_size = 1000;

  redshift_range[0] = 0.0;
  redshift_range[1] = 10.0;

  if(lf != NULL)
    luminosity_function = lf;
  else
    cout << "ERROR: NULL Pointer Passed to Simulator" << endl;
  
  models = new model_lib(modfile);
  observations = new obs_lib(obsfile);
  
  for (int i=0;i<3;i++){
    bands[i] = b[i];
    band_errs[i] = b_err[i];
    flux_limits[i]= f_lims[i];
  }
}

void simulator::set_bands(double b[],double b_err[],double f_lims[]){
  for (int i=0;i<3;i++){
    bands[i] = b[i];
    band_errs[i] = b_err[i];
    flux_limits[i]= f_lims[i];
  }
}

void simulator::set_lumfunct(lumfunct *lf){
  if(lf != NULL)
    luminosity_function = lf;
  else
    cout << "ERROR: NULL Pointer Passed to Simulator" << endl;
}

void simulator::set_model_lib(string modfile){
  models = new model_lib(modfile);
}

void simulator::set_obs_lib(string obsfile){
  observations = new obs_lib(obsfile);
}

void simulator::set_zrange(double zrange[]){
  if(zrange[0] > 0.0)
    redshift_range[0] = zrange[0];
  if(zrange[1] < 10.0)
    redshift_range[1] = zrange[1];
}

void simulator::set_simulated_size(double size){
  distribution_size = size;
}

double simulator::simulate(const gsl_vector *v,void *params){
  if(models != NULL){
    static fixed_params *p;
    p = (fixed_params *)params;
    static int pnum(p->pnum+1);

    double mean[pnum];
    double sigma[pnum];
    double min[pnum];
    double max[pnum];
    static double alt_size = 0;

    mean[0] = gsl_vector_get(v,0);
    sigma[0] = abs(gsl_vector_get(v,1));
    min[0] = redshift_range[0];
    max[0] = redshift_range[1];
    
    for (int i=1;i<=p->pnum;i++){
      mean[i] = p->p[i-1].mean;
      sigma[i] = p->p[i-1].sigma;
      min[i] = p->p[i-1].min;
      max[i] = p->p[i-1].max;
    }
    
    if(distributions == NULL)
      distributions = new sim(models);
    distributions->initialize_dists(mean,sigma,distribution_size,min,max);

    if(p->znum > 1){
      alt_size = distribution_size*abs(gsl_vector_get(v,2));
      mean[0] = gsl_vector_get(v,3);
      sigma[0] = abs(gsl_vector_get(v,4));
      if(alt_size > 0)
	distributions->add_dists(mean,sigma,alt_size,min,max);
    }

    if(bands[0] != -1){
      distributions->simulate(bands,band_errs,flux_limits,*luminosity_function);
    }
    else{
      cout << "ERROR: Bands Uninitialized" << endl;
      return -1;
    }
    
    double *c1,*c2,*mc1,*mc2;
    int osize = observations->get_snum();
    int msize = distributions->size();
    if (diagnostic == NULL)
      diagnostic = new hist_lib;

    observations->get_all_colors(c1,c2);
    distributions->get_all_colors(mc1,mc2);

    diagnostic->init_obs(c1,c2,osize);
    diagnostic->init_model(mc1,mc2,msize);
    
    delete[] c1;
    delete[] c2;
    delete[] mc1;
    delete[] mc2;

    return diagnostic->get_chisq();
  }
  else{
    cout << "ERROR: NULL Model Library" << endl;
    return -1;
  }
}

bool simulator::save(string outfile){
  bool hist = diagnostic->write_fits(outfile);
  
  using namespace CCfits;
  
  std::auto_ptr<FITS> pfits(0);

  try{
    pfits.reset(new FITS(outfile,Write));
  }
  catch (CCfits::FITS::CantOpen){
    return false;       
  }
  
  int pnum = models->get_pnum();
  
  int size = distributions->size();
  valarray<double> f1(size),f2(size),f3(size),luminosity(size);
  valarray<double> redshift(size);
  valarray<double> param(size);
  valarray<double> params[pnum];
  for (int i=0;i<pnum;i++)
    params[i].resize(size);
  static sprop temp;
  double ptemp[pnum];
  static double extra_z;

  for(int i=0;i<size;i++){
    temp = distributions->get(i);
    distributions->get(i,extra_z,ptemp);
    f1[i] = temp.fluxes[0];
    f2[i] = temp.fluxes[1];
    f3[i] = temp.fluxes[2];
    param[i] = temp.mod_num;
    redshift[i] = temp.redshift;
    luminosity[i] = temp.luminosity;
    for(int j=0;j<pnum;j++)
      params[j][i] = ptemp[j];
  }

  static std::vector<string> colname((6+pnum),"");
  static std::vector<string> colunit((6+pnum),"");
  static std::vector<string> colform((6+pnum),"f13.8");
  
  static string pnames[] = {"P0","P1","P2","P3","P4"};

  colname[0] = "F1";
  colname[1] = "F2";
  colname[2] = "F3";
  colname[3] = "Z";
  colname[4] = "M";
  colname[5] = "Lum";
  
  colunit[0] = "Jy";
  colunit[1] = "Jy";
  colunit[2] = "Jy";
  colunit[3] = "-";
  colunit[4] = "-";
  colunit[5] = "W/Hz";

  colform[5] = "e13.5";

  for(int i=6;i<6+pnum;i++){
    colname[i] = pnames[i-6];
    colunit[i] = "-";
  }

  static string hname("Parameter Distributions");
  Table *newTable= pfits->addTable(hname,size,colname,colform,colunit,AsciiTbl);

  newTable->column(colname[0]).write(f1,1);
  newTable->column(colname[1]).write(f2,1);
  newTable->column(colname[2]).write(f3,1);
  newTable->column(colname[3]).write(redshift,1);
  newTable->column(colname[4]).write(param,1);
  newTable->column(colname[5]).write(luminosity,1);
  for (int i=6;i<6+pnum;i++)
    newTable->column(colname[i]).write(params[i-6],1);

  return hist;
}

simulator::~simulator(){

  if(distributions != NULL)
    delete distributions;

  if(models != NULL)
    delete models;
  
  if(observations != NULL)
    delete observations;

  if(diagnostic != NULL)
    delete diagnostic;
}

