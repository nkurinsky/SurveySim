#include "hist_lib.h" //See this header file for additional function description

hist_lib::hist_lib(){
  model_hist = NULL;
  obs_hist = NULL;
  comparison_hist = NULL;
  range[0] = -1;
  range[1] = -1;
  binsize = -1;
  xysize = -1;
}

hist_lib::hist_lib(double obs_c1[],double obs_c2[],int obs_size){
  static double temp_range[2];
  static double x,y;
 
  //set the range to the absolute extremes of the color arrays
  x = gsl_stats_min(obs_c1,H_STRIDE,obs_size);
  y = gsl_stats_min(obs_c2,H_STRIDE,obs_size);
  if (x <= y)
    temp_range[0] = x;
  else
    temp_range[0] = y;

  x = gsl_stats_max(obs_c1,H_STRIDE,obs_size);
  y = gsl_stats_max(obs_c2,H_STRIDE,obs_size);
  if (x >= y)
    temp_range[1] = x;
  else
    temp_range[1] = y;

  //create histogram for observation and store it
  obs_hist = get_hist(obs_c1,obs_c2,obs_size,temp_range);

  //initialize to null to indicate that it has yet to be constructed
  model_hist = NULL;
  comparison_hist = NULL;
  osize = obs_size;
  msize = -1;
}

void hist_lib::init_obs(double c1[],double c2[],int size){
  if(obs_hist != NULL){
    for (int i=0;i<xysize;i++){
      delete[] obs_hist[i];
    }
    delete[] obs_hist;
    obs_hist = NULL;
  }
  static double temp_range[2];
  static double x,y;
  
  //set the range to the absolute extremes of the color arrays
  x = gsl_stats_min(c1,H_STRIDE,size);
  y = gsl_stats_min(c2,H_STRIDE,size);
  if (x <= y)
    temp_range[0] = x;
  else
    temp_range[0] = y;
  
  x = gsl_stats_max(c1,H_STRIDE,size);
  y = gsl_stats_max(c2,H_STRIDE,size);
  if (x >= y)
    temp_range[1] = x;
  else
    temp_range[1] = y;
  
  //create and store histogram
  obs_hist = get_hist(c1,c2,size,temp_range);
  osize = size;

  if(comparison_hist != NULL){
    for (int i=0;i<xysize;i++){
      delete[] comparison_hist[i];
    }
    delete[] comparison_hist;
    comparison_hist = NULL;
  }
}

bool hist_lib::init_model(double c1[],double c2[],double weights[],int size){
  if (model_hist != NULL){
    for (int i=0;i<xysize;i++){
      delete[] model_hist[i];
    }
    delete[] model_hist;
    model_hist = NULL;
  } 
  //create histogram with class parameters
  if(obs_hist != NULL){
    model_hist = compute_hist(c1,c2,weights,size);
    msize = size;
    
    if(comparison_hist != NULL){
      for (int i=0;i<xysize;i++){
	delete[] comparison_hist[i];
      }
      delete[] comparison_hist;
      comparison_hist = NULL;
    }
  }
  return true;
}

double hist_lib::get_chisq(){
  if(comparison_hist != NULL)
    return chisq;
  else if ((obs_hist != NULL) and (model_hist != NULL)){
    chisq = fit_err();
    return chisq;
  }
  
  printf("%s \n","ERROR: (Hist) Model and Observation not fully initialized");
  return -1;
}

bool hist_lib::write_fits(string filename){
  filename = "!"+filename;
  
  static long naxis(2);
  static long naxes[] = {xysize,xysize};

  using namespace CCfits;
  std::auto_ptr<FITS> pFits(0);
  try{
    pFits.reset(new FITS(filename,DOUBLE_IMG,naxis,naxes));
  }
  catch (FITS::CantCreate){
    return false;
  }
  
  static std::vector<long> extAx(2,xysize);
  static string names[] = {"model","observation"};
  static int nelements;
  nelements = xysize*xysize;

  static std::valarray<double> model_1d(nelements);
  static std::valarray<double> obs_1d(nelements);
  static std::valarray<double> comparison_1d(nelements);

  for (int i=0;i<xysize;i++){
    for (int j=0;j<xysize;j++){
      model_1d[i*xysize+j] = model_hist[i][j];
      obs_1d[i*xysize+j] = obs_hist[i][j];
      comparison_1d[i*xysize+j] = comparison_hist[i][j];
    }
  }

  pFits->pHDU().addKey("DIM",xysize,"Side length of square histogram area"); 
  pFits->pHDU().addKey("H_MIN",range[0],"Minimum Histogram Bin Value");
  pFits->pHDU().addKey("H_MAX",range[1],"Maximum Histogram Bin Value");
  pFits->pHDU().addKey("BINSIZE",binsize,"Histogram bin size");
  pFits->pHDU().addKey("FITSTAT",chisq,"Fitting Statistic");

  pFits->pHDU().write(1,nelements,comparison_1d);

  ExtHDU* modImage = pFits->addImage(names[0],LONG_IMG,extAx);
  ExtHDU* obsImage = pFits->addImage(names[1],LONG_IMG,extAx);

  modImage->write(1,nelements,model_1d);
  obsImage->write(1,nelements,obs_1d);

  return true;
}

hist_lib::~hist_lib(){
  if(model_hist != NULL){
    for (int i=0;i<xysize;i++){
      delete[] model_hist[i];
    }
    delete[] model_hist;
  }
  if(obs_hist != NULL){
    for (int i=0;i<xysize;i++){
      delete[] obs_hist[i];
    }
    delete[] obs_hist;
  }
  if(comparison_hist != NULL){
    for (int i=0;i<xysize;i++){
      delete[] comparison_hist[i];
    }
    delete[] comparison_hist;
  }
}

double ** hist_lib::get_hist(double c1[],double c2[],int cnum,double inp_range[]){
  static double ** temp;
  static double sd1,sd2;
  sd1 = gsl_stats_sd(c1,H_STRIDE,cnum);
  sd2 = gsl_stats_sd(c2,H_STRIDE,cnum);

  //compute binsize from Scott 1979 formalism (binsize=3.49*sd*n^(-1/3))
  if (sd1 >= sd2)
    binsize = 3.49*sd1/pow(cnum,0.33333);
  else
    binsize = 3.49*sd2/pow(cnum,0.33333);

  range[0] = inp_range[0];
  range[1] = inp_range[1];
  static double rdiff; 
  rdiff = (range[1]-range[0]);
  xysize = ceil(rdiff/binsize); //ensures that range is whole multiple of binsize
  range[1] = range[0]+xysize*binsize; //recompute upper bound
  
  double weights[cnum];
  for (int i=0;i<cnum;i++)
    weights[i] = 1.0;

  //create histogram with given range and calculated binsize
  temp = compute_hist(c1,c2,weights,cnum);
  return temp;
}

double ** hist_lib::compute_hist(double c1[],double c2[],double weights[],int cnum){
  static double ** retvals;

  //initialize two dimensional integer array (Dynamic)
  retvals = new double*[xysize];
  for (int i=0;i<xysize;i++){
    retvals[i] = new double[xysize];
    for (int j=0;j<xysize;j++){
      retvals[i][j] = 0;
    }
  }
  
  //computes indices of point using integer rounding. increments integer at associated
  //index (eq used: (point-min)/binsize   round down to get index)
  static int xpt,ypt;
  static double exc;
  xpt = 0;
  ypt = 0;
  exc = 0;
  for (int i=0;i<cnum;i++){
    if((c1[i] >= range[0]) and (c1[i] < range[1]) and (c2[i] >= range[0]) and (c2[i] < range[1])){
      xpt = int((c1[i]-range[0])/binsize);
      ypt = int((c2[i]-range[0])/binsize);
      retvals[xpt][ypt]+= weights[i];
    }
    else
      exc+= weights[i];
  }
  m_exc = exc;
  //number of points which lay outside of the specified range
  printf("CHist Exclusions: %6.2E %6.2E (%5.2f%%)\n",double(exc),double(cnum),(double(exc)/double(cnum))*100.0);
  
  return retvals;
}

//here replace sqrt with tabulated vals and new formula
double hist_lib::fit_err(){
  //Chi-square
  static double chisq,obs_temp,mod_temp,obs_temp_err,mod_temp_err,err_temp,chisq_temp;
  static bool obs_valid,mod_valid;
  
  comparison_hist = new double*[xysize];
  for (int i=0;i<xysize;i++){
    comparison_hist[i] = new double[xysize];
  }

  chisq = 0;

  for (int i=0;i<xysize;i++){
    for (int j=0;j<xysize;j++){
      obs_temp = double(obs_hist[i][j]);
      mod_temp = double(model_hist[i][j]);

      obs_valid = (obs_temp > 0);
      mod_valid = (mod_temp > 0);

      chisq_temp = 0;
      if(obs_valid or mod_valid){
	if(mod_valid){
	  //mod_temp_err = sqrt(mod_temp);
	  mod_temp_err = poiss_err(mod_temp);
	  mod_temp /= double(msize);
	  mod_temp_err /= double(msize);
	}
	else{
	  mod_temp = 0;
	  mod_temp_err = 1.841;
	}

	if(obs_valid){
	  obs_temp_err = poiss_err(obs_temp);
	  obs_temp /= double(osize);
	  obs_temp_err /= double(osize);
	}
	else{
	  obs_temp = 0;
	  obs_temp_err = 1.841;
	}
       
	err_temp = sqrt(pow(mod_temp_err,2)+pow(obs_temp_err,2));
	chisq_temp = (pow((obs_temp-mod_temp),2)/pow(err_temp,2));
	chisq += chisq_temp;
      }	
      //comparison_hist[i][j] = int((double(model_hist[i][j])*double(osize))/double(msize))-obs_hist[i][j];
      comparison_hist[i][j] = chisq_temp;
    }
  }
  return chisq;
}

double hist_lib::poiss_err(int x){
  if (x le 0)
    return poiss[0];
  else if (x le 100)
    return poiss[x];
  else
    return sqrt(x)+1.0;
}
