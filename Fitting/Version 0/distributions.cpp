#include "distributions.h" //See this header for more information on the following funcitons

distribution::distribution(gsl_rng * r,int size,double mean,double sigma,double range[]){
  double * result = gauss_random(r,range,mean,sigma,size);
  values.assign(result,result+size);
  delete[] result;
}

distribution::distribution(gsl_rng * r,int size,double range[]){
  double * result = random(r,range,size);
  values.assign(result,result+size);
  delete[] result;
}

distribution::distribution(int size,double value){
  for (int i=0;i<size;i++)
    values.push_back(value);
}

void distribution::add(gsl_rng * r,int size,double mean,double sigma,double range[]){
  double * result = gauss_random(r,range,mean,sigma,size);
  for (int i=0;i<size;i++){
    values.push_back(result[i]);
  }
  delete[] result;
}

void distribution::add(gsl_rng * r,int size,double range[]){
  double * result = random(r,range,size);
  for (int i=0;i<size;i++){
    values.push_back(result[i]);
  }
  delete[] result;
}

void distribution::add(int size,double value){
  for (int i=0;i<size;i++)
    values.push_back(value);
}

double distribution::get(int i){
  if(i < int(values.size()))
    return values[i];
  else{
    printf("%s \n","ERROR: (Dists) Index Out of Redshift Array Bounds");
    return -1.0;
  }
}

void distribution::get_all(double retvals[]){
  for (int i=0;i<int(values.size());i++){
    retvals[i] = values[i];
  }
  return;
}

double lf(double x,void * params){
  lf_params *p = (lf_params *)params;

  double term1 = p->OmegaM*pow((1+x),3);
  double result = 1/sqrt(term1+p->Lambda0);
  return result;
}

lumfunct::lumfunct(double flimit){
  flux_limit = flimit;
  phi0 = pow(10,-2.2);
  L0 = 23.69;
  alpha = 0.47;
  beta = 2.88;
  p = -6.7;
  q = 3.5;
  exp_range[0] = -3;
  exp_range[1] = 5;
  zmax = 10.0000e0;
  H0 = 73;
  Lambda0 = 0.728;
  OmegaM = 0.272;
  init = false;

  gsl_rng_default_seed = time(NULL);
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
}

double lumfunct::lumdist(double z){
  return gsl_spline_eval(ldist_spline,z,ldist_acc);
}

void lumfunct::set_phi0(double val){
  phi0=val;
}

void lumfunct::set_L0(double val){
  L0=val;
}

void lumfunct::set_alpha(double val){
  alpha=val;
}

void lumfunct::set_beta(double val){
  beta=val;
}

void lumfunct::set_p(double val){
  p=val;
}

void lumfunct::set_q(double val){
  q=val;
}

void lumfunct::set_flux_limit(double val){
  if(val > 0)
    flux_limit=val;
  else
    printf("%s \n","ERROR: Invalid Flux Limit");
}

void lumfunct::set_range(double lower,double upper){
  if(lower < upper){
    exp_range[0] = lower;
    exp_range[1] = upper;
  }
  else
    printf("%s \n","ERROR: Range Invalid");
}

void lumfunct::set_H0(double val){
  if(val > 0)
    H0 = val;
  else
    printf("%s \n","ERROR: Invalid H0 Value");
}

void lumfunct::set_Lambda0(double val){
  if (val > 0)
    Lambda0 = val;
  else
    printf("%s \n","ERROR: Invalid Lamda0 Value");
}

void lumfunct::set_OmegaM(double val){
  if (val > 0)
    OmegaM = val;
  else
    printf("%s \n","ERROR: Invalid OmegaM Value");
}

double lumfunct::get_amp(double redshift, double &lum){
  double dist = lumdist(redshift);
  int ind = int(redshift*10);
  double flux;
  int cnum = 0;
  do{
    lum = luminosities[ind][gsl_ran_discrete(r,g_array[ind])];
    flux = (lum*1e26)/(4*M_PI*pow(dist,2)*(1+redshift));
    cnum++;
  }while (flux < flux_limit);
  
  return flux;
}

void lumfunct::initialize(){

  if(init){
    delete [] z_ldists;
    delete [] ldists;
    gsl_spline_free(ldist_spline);
    gsl_interp_accel_free(ldist_acc);
  }

  gsl_integration_workspace * w;
  double result, error;
  
  lf_params params = {OmegaM,Lambda0};
  
  gsl_function F;
  F.function = &lf;
  F.params = &params;
  
  double omegaK = 1.0-Lambda0-OmegaM;
  int arraysize = zmax*1000+1;
  ldists = new double[arraysize];
  z_ldists = new double[arraysize];
  
  z_ldists[0] = 0.0;
  ldists[0] = 0.0;
  int mi = 1;
  
  if(abs(omegaK) < 0.01){
    for(double z=1e-3;z<=zmax;z+=1e-3){
      result = 0;
      error = 0;
      w = gsl_integration_workspace_alloc(1000);
      gsl_integration_qags (&F,(z-(1e-3)),z,0,1e-7,1000,w,&result,&error); 
      gsl_integration_workspace_free (w);
      
      result *= (GSL_CONST_MKSA_SPEED_OF_LIGHT*(1+z)/(1e3*H0))*MPC_TO_METER;
      z_ldists[mi] = z;
      ldists[mi] = ldists[mi-1] + result;
      mi++;
    }
  }
  else{
    printf("%s \n","ERROR: Lumdist cannot handle curvature");
    result = -1;
  }
  
  ldist_acc = gsl_interp_accel_alloc();
  ldist_spline = gsl_spline_alloc(gsl_interp_cspline,arraysize);

  gsl_spline_init(ldist_spline,z_ldists,ldists,arraysize);

  int izmax = int(zmax*10.0);

  //deallocate memory if already allocated so that it make be reallocated
  if(init){
    for (int i=0;i<izmax;i++){
      gsl_ran_discrete_free(g_array[i]);
      delete[] luminosities[i];
    }
    delete[] g_array;
    delete [] luminosities;
  }

  //allocate memory
  int lumsize = 1000;
  double ratios[lumsize];
  luminosities = new double*[izmax];
  for (int i=0;i<izmax;i++)
    luminosities[i] = new double[lumsize];

  //initialize variables for function computation
  double lumfunct[lumsize];
  double step_size = (exp_range[1]-exp_range[0])/double(lumsize-1);

  double ratio,scale,dist,lmin_exp;
  int ind,zind;

  g_array = new gsl_ran_discrete_t*[izmax];
  
  zind = 0;
  for (double z=0.000;z<(zmax-0.05);z+= 1e-1){
    ind = 0;
    dist = lumdist(z); //distance for given redshift
    //find mimimum value for L/L* from redshift and flux limit, with K correction
    lmin_exp = log10(flux_limit*1e-26*4*M_PI*pow(dist,2)*(1+z))-L0-log10(pow((1+z),q));
    for (double i=exp_range[0];i<=exp_range[1];i+=step_size){
      ratio = pow(10,i);
      scale = phi0*pow((1+z),p);
      ratios[ind] = ratio;
      luminosities[zind][ind] = ratio*pow(10,L0)*pow((1+z),q);
      //If ratio is lower than lower bound, assign it a 0 detection probability, otherwise compute function
      if (i > lmin_exp)
	lumfunct[ind] = scale/(pow(ratio,alpha)+pow(ratio,beta));
      else
	lumfunct[ind] = 0.00;
      ind++;
    }
    //initialize random number lookup table
    g_array[zind] = gsl_ran_discrete_preproc(lumsize,lumfunct);
    zind++;
  }

  init = true;
}

lumfunct::~lumfunct(){

  for (int i=0;i<int(zmax*10);i++){
    gsl_ran_discrete_free(g_array[i]);
    delete[] luminosities[i];
  }
  delete[] luminosities;
  delete[] g_array;
  gsl_rng_free(r);

  delete [] z_ldists;
  delete [] ldists;
  
  gsl_spline_free(ldist_spline);
  gsl_interp_accel_free(ldist_acc);
}
