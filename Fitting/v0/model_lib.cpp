#include "model_lib.h" //see header for function descriptions and additional comments

model::model(double *bands,double * fluxes,int fluxnum){
  this->fluxes.assign(fluxes,fluxes+fluxnum);
  
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline,fluxnum);
  gsl_spline_init(spline,bands,fluxes,fluxnum);
}

double model::get_flux(double band,double redshift){
  double new_band = band/(1+redshift);
  return gsl_spline_eval(spline,new_band,acc);
}

model::~model(){
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}

model_lib::model_lib(string fitsfile){
  //initialize FITS input
  FITS* pInfile;
  pInfile = new FITS(fitsfile,Read);
  
  //reading primary header from fits file
  HDU& table = pInfile->pHDU();
  //extracting conents of image from fits file
  std::valarray<double> contents;
  PHDU& image = pInfile->pHDU();
  image.read(contents);

  string ti[10] = {"0","1","2","3","4","5","6","7","8","9"};
  double scale = 1;
  double temp;

  //get number of parameters and initialize dynamic memory
  table.readKey("P_NUM",pnum);
  psize = new int[pnum];

  //get min, max, and step size for each parameter and compute
  //total number of parameter values used in model library
  for (int i=0;i<pnum;i++){
    table.readKey("P"+ti[i]+"_SIZE",temp);
    psize[i] = int(temp);
    scale *= temp; //gives number of models total
  }

  //generate wavelength vector from fits keywords
  double wave_min,wave_max,wave_step;
  table.readKey("WAVE_MIN",wave_min);
  table.readKey("WAVE_MAX",wave_max);
  table.readKey("WAVE_SEP",wave_step);
  band_num = ((wave_max-wave_min)/wave_step)+1;

  bands = new double[band_num];

  int ind = 0;
  for (double i=wave_min;i <= (wave_max+(0.5*wave_step));i+=wave_step){
    bands[ind] = i;
    ind++;
  }
  printf("\nModel Wavelength Range: %9.3e - %9.3e\nNumber of Data Points: %i\n",bands[0],bands[band_num-1],band_num);
  
  int n_fluxes = image.axis(0)/scale;
  double * fluxes = new double[n_fluxes];
  model * new_model;
  //int * ids;
  //double * params;

  mod_num = scale;
  models = new model*[int(scale)];

  //initialize models by extracting fluxes from overall table
  for (int mod=0;mod < scale;mod++){
    for (int fi=0;fi<n_fluxes;fi++){
      fluxes[fi] = contents[fi*scale+mod];
    }
    //ids = get_ids(mod);
    //params = get_params(ids);
    new_model = new model(bands,fluxes,n_fluxes);
    models[mod]=new_model;
    //delete[] ids;
    //delete[] params;
  }

  acc = gsl_interp_accel_alloc ();
  interp = gsl_interp_alloc (gsl_interp_cspline, 3);

  //delete unwanted dynamic variables
  delete[] fluxes;
  delete pInfile;
}

double model_lib::get_flux(double ids[],double band,double redshift){
  static double modnums[3];
  int ptemp[pnum];
  static double fluxes[3];
  static int up,down;

  for (int i=0;i<3;i++){
    for (int j=0;j<pnum;j++){
      up = ceil(ids[j]);
      down = floor(ids[j]);
      if(up == down)
	ptemp[j] = up;
      else if(down == 0)
	ptemp[j] = i;
      else if(up == (psize[j]-1))
	ptemp[j] = psize[j]-3+i;
      else
	ptemp[j] = down-1+i;
    }
    modnums[i] = index(ptemp);
    fluxes[i] = models[int(modnums[i])]->get_flux(band,redshift);
  }
  
  gsl_interp_accel_reset(acc);
  gsl_interp_init (interp, modnums, fluxes, 3);
  return gsl_interp_eval (interp,modnums,fluxes,index(ids),acc);
}

int model_lib::index(int params[]){
  int retval = 0;
  int scale = 1;
  for (int i=0;i<pnum;i++){
    if (i > 0) 
      scale*=psize[i-1];
    retval += params[i]*scale;
  }
  return retval;
}

double model_lib::index(double params[]){
  double retval = 0;
  double scale = 1;
  for (int i=0;i<pnum;i++){
    if (i > 0) 
      scale*=psize[i-1];
    retval += params[i]*scale;
  }
  return retval;
}

int * model_lib::get_ids(int index){
  int * retvals = new int[pnum];
  for (int i=0;i<pnum;i++){
    retvals[i] = index % psize[i];
    index = index/psize[i];
  }
  return retvals;
}

void model_lib::get_psize(int ps[]){
  for (int i=0;i<pnum;i++){
    ps[i] = psize[i];
  }
}

model_lib::~model_lib(){
  delete[] psize;
  for (int i=0;i<mod_num;i++)
    delete models[i];
  delete[] models;
  delete[] bands;

  gsl_interp_free (interp); 
  gsl_interp_accel_free (acc);
}
