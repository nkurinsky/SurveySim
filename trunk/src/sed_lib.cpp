#include "sed_lib.h"

sed::sed(){
  interp_init = false;
  zdim=false;
}

sed::sed(double * f,double *bands, int bnum) : sed(){
  
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline,bnum);
  gsl_spline_init(spline,bands,f,bnum);
  
  brange[0] = bands[0];
  brange[1] = bands[bnum-1];
  
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
      printf("ERROR: emitted lambda %e out of range %e -> %e\n",x_em,brange[0],brange[1]);
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

sed_lib::sed_lib(string fitsfile, int nz, double zmin, double dz){

  redshiftq=false;
  typeq=false;
  interp_init=false;

  color_exp = 0;
  color_zcut = 0;
  color_evolution = 1;

  interp_znum=nz;
  interp_zmin=zmin;
  interp_dz=dz;
  
  unique_ptr<FITS> pfits;
  try{
    pfits.reset(new FITS(fitsfile));
  }
  catch(CCfits::FITS::CantOpen){
    cout << "Cannot open "<< fitsfile << endl;
    exit(1);
  }

  HDU& info=pfits->pHDU();
  info.readAllKeys();
  auto settings=info.keyWord();

  int nmodels,nzbins;
  if(settings.count("NTYPES") == 1){
    settings["NTYPES"]->value(nmodels);
  }
  else{
    cout << "Keyword NTYPES not specified; exiting" << endl;
    exit(1);
  }

  if(settings.count("NZBINS") == 1){
    settings["NZBINS"]->value(nzbins);
  }
  else{
    cout << "Keyword NZBINS not specified; exiting" << endl;
    exit(1);
  }

  float scale;
  if(settings.count("SCALE") == 1){
    settings["SCALE"]->value(scale);
  }
  else{
    cout << "Keyword SCALE not specified; exiting" << endl;
    exit(1);
  }

  redshiftq=nzbins>1;
  typeq=nmodels>1;

  vector<double> zmins;
  vector<double> zmaxs;
  vector<string> types;

  //read SED types
  string type;
  for(int i=1;i<=nmodels;i++){
    char buffer[10];
    sprintf(buffer,"SEDTYPE%i",i);
    string field(buffer);
    if(settings.count(field) > 0){
      settings[field]->value(type);
      types.push_back(type);
    }
    else{
      cout << "Keyword " << field << " not specified; exiting" << endl;
      exit(1);
    }
  }

  //Read redshift bins
  double tzmin;
  for(int i=0;i<nzbins;i++){
    char buffer[10];
    sprintf(buffer,"ZMIN%i",i);
    string field(buffer);
    if(settings.count(field) > 0){
      settings[field]->value(tzmin);
    }
    else{
      printf("Keyword %s missing from header (%s)\n",field.c_str(),fitsfile.c_str());
      exit(1);
    }
    zmins.push_back(tzmin);    
  }

  //Read redshift bins
  double tzmax;
  for(int i=0;i<nzbins;i++){
    char buffer[10];
    sprintf(buffer,"ZMAX%i",i);
    string field(buffer);
    if(settings.count(field) > 0){
      settings[field]->value(tzmax);
    }
    else{
      printf("Keyword %s missing from header (%s)\n",field.c_str(),fitsfile.c_str());
      exit(1);
    }
    zmaxs.push_back(tzmax);
  }

  auto extensions=pfits->extension();

  //read luminosities from primary header
  lnum=extensions.begin()->second->numCols()-1;
  lums.reset(new double[lnum]);
  unique_ptr<double[]> inds(new double[lnum]);
  char buffer[10];
  double temp;
  for (unsigned int i=0;i<lnum;i++) {
    std::string a = "L";
    sprintf(buffer,"%i",i+1);
    a.append(buffer);
    if(settings.count(a) == 1){
      settings[a]->value(temp);
      lums[i] = temp;
      inds[i] = static_cast<double>(i);
    }
    else{
      printf("Keyword %s missing from header (%s)\n",a.c_str(),fitsfile.c_str());
      exit(1);
    }
  }
  
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline,lnum);
  gsl_spline_init(spline,lums.get(),inds.get(),lnum);

  //find extensions corresponding to types
  map<string,vector<string> > type_exts;
  for(auto itype=types.begin();itype!=types.end();++itype){
    for (auto ext=extensions.begin();ext != extensions.end();++ext){
      std::size_t found = ext->first.find(*itype);
      if (found!=std::string::npos){
	type_exts[*itype].push_back(ext->first);
      }      
    }
  }
  tnum=type_exts.size();
  
  cout << "Finding Types" << endl;
  unique_ptr<sed> new_sed;  
  LOG_INFO(printf("Reading in Types:\n"));
  //for (auto itr=type_exts.begin(); itr != type_exts.end(); ++itr){
  for(auto itr=types.begin();itr != types.end();++itr){
    if(type_exts.count(*itr) == 0)
      continue;

    int type_len=itr->length();
    vector<double> lambda;
    vector<double> zs;
    for(auto eitr=type_exts[*itr].begin();eitr!=type_exts[*itr].end();++eitr){
      string extstr(*eitr);
      extstr.erase(0,type_len+2);
      int zbin(atoi(extstr.c_str()));
      zs.push_back(zmins[zbin]);
    }

    auto firstInstance=extensions.find(type_exts[*itr].front())->second;
    int tablelength=firstInstance->rows();
    firstInstance->column(1).read(lambda,0,tablelength);
    static double scale_factor(pow(10,scale));
    for(auto lam=lambda.begin();lam!=lambda.end();++lam)
      *lam *= scale_factor;
    LOG_INFO(printf("%s\t",(*itr).c_str()));
    if(itr == types.begin()){
      brange[0] = lambda.front();
      brange[1] = lambda.back();
    }
    LOG_INFO(printf("\n"));
      
    for(int i=2;i<lnum+2;i++){
      vector<double> fluxtemp;
      vector<double> fluxes;
      for(auto eitr=type_exts[*itr].begin();eitr!=type_exts[*itr].end();eitr++){
	extensions.find(*eitr)->second->column(i).read(fluxtemp,0,tablelength);
	for(int i=0;i<fluxtemp.size();i++){
	  fluxes.push_back(fluxtemp[i]);
	}
      }
      if(zs.size() > 1)
	new_sed.reset(new sed(fluxes.data(),
			      lambda.data(),
			      lambda.size(),
			      zs.data(),
			      zs.size()));
      else
	new_sed.reset(new sed(fluxes.data(),
			      lambda.data(),
			      lambda.size()));
      seds.push_back(move(new_sed));
    }
  }

  color_exp = 0;
}

bool sed_lib::load_filters(string file,int lflag_tmp){
  
  if(filters.load_filters(file,lflag_tmp)){
    logflag=lflag_tmp;
    w = gsl_integration_workspace_alloc(SL_INT_SIZE);
    return true;
  }
  
  return false;
}

void sed_lib::get_filter_info(string names[], double limits[], double errors[]){
  if (filters.init()){
    filters.filter_info(names,limits,errors);
  }
  else{
    printf("ERROR: Filter library not initialized!!!\n");
    exit(1);
  }
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
      return seds[i]->get_flux(band,redshift);
    else{
      LOG_CRITICAL(printf("%s\n%s\n","ERROR: Band out of model range.","Check that the obs bands are within range and unit consistent"));
    }
  }
  else{
    LOG_CRITICAL(printf("%s %i \n","ERROR: Unitialized Model Called!!!",i));
  }
  
  return -1;
}

double sed_lib::get_filter_flux(double lum, double redshift, short sedtype, short filter_id){

  if((sedtype < 0) or (sedtype >= tnum)){
    printf("ERROR: Invalid sedtype \"%i\" out of possible range (0,%i)\n",sedtype,tnum);
    return -1;
  }
  else if((filter_id < 0) or (filter_id > 2)){
    printf("ERROR: Filter number %i invalid, 0-2 acceptable.\n",filter_id);  
    return -1;
  }
  else if (not filters.init()){
    printf("ERROR: Filter library not initialized!!!\n");
    return -1;
  }
   
  if(not interp_init)
    initialize_filter_fluxes(logflag);
  
  return interpolate_flux(lum,redshift,sedtype,filter_id);
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
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  
  if(filters.init())
    gsl_integration_workspace_free(w);
}

void sed_lib::initialize_filter_fluxes(int logflag){
  //vector contains alg_lib interpolation class
  flux_interpolator.resize(FILTER_NUM*tnum);
  
  alglib::real_1d_array ls;
  alglib::real_1d_array zs;
  alglib::real_1d_array fs;

  ls.setcontent(lnum,lums.get());
  zs.setlength(interp_znum);
  for(int zi=0;zi<interp_znum;zi++)
    zs[zi] = static_cast<double>(zi)*interp_dz+interp_zmin;
  fs.setlength(lnum*interp_znum);
  
  for (int type=0;type<tnum;type++){
    for (int filter=0;filter<3;filter++){
      for(int zi=0;zi<interp_znum;zi++){  
	for(int li=0;li<lnum;li++){
	  fs[li+zi*lnum]=convolve_filter(li,zs[zi],type,filter);
	}
      }
      LOG_DEBUG(printf("\rBuilding Spline..."));
      try{
	alglib::spline2dbuildbilinearv(ls,lnum,
				       zs,interp_znum,
				       fs,1,
				       flux_interpolator[type*FILTER_NUM+filter]);
      }
      catch(alglib::ap_error e){
	printf("\nERROR: Failed to create interpolation for Filter Convolution\n");
	printf("error msg: %s\n", e.msg.c_str());
	exit(1);
      }  
    }
  }
  
  LOG_INFO(printf("\rBuilt Spline Interpolation for Filter Fluxes\n"));
  interp_init=true;
}

double sed_lib::interpolate_flux(double lum, double redshift, short sedtype, short filter_id){
  double retval(-1);
  try{
    retval =  alglib::spline2dcalc(flux_interpolator[sedtype*FILTER_NUM+filter_id],lum,redshift);
    
    if(redshift < color_zcut)
      retval *= pow((1.0+redshift),color_exp);
    else
      retval *= color_evolution;
  }
  catch(alglib::ap_error e){
    printf("ERROR: Failed to interpolate for z=%lf, lum=%lf, sed=%i, filter=%i\n",redshift,lum,sedtype,filter_id);
    printf("error msg: %s\n", e.msg.c_str());
  }
  return retval;
}

double sed_lib::convolve_filter(short lum_id, double redshift, short sedtype, short filter_id){

  static map<tuple<short,short,double>,double > fluxes[3];
  tuple<short,short,double> params (lum_id,sedtype,redshift);

  if(filters.get(filter_id).low() >= filters.get(filter_id).high()){
    LOG_CRITICAL(printf("ERROR: Filter %i unititialized.\n",filter_id));
    return -1;
  }

  if((filters.get(filter_id).low() < brange[0]) or (filters.get(filter_id).high() > brange[1])){
    LOG_CRITICAL(printf("ERROR: Filter out of model range.\nCheck that the obs bands are within range and unit consistent\n"));
    return -1;
  }

  //we use a map here to ensure that a given redshift/SED/filter combination is only convolved once; we recompute color evolution every time as this may change in fitting, while filters/SEDs are static.
  if(fluxes[filter_id].count(params) == 0){
    double scale,result,error;
    double bounds[2];
    gsl_function F;
    flux_yield_params p(seds[lum_id+lnum*sedtype].get(),&filters.get(filter_id),redshift);
    F.function = &flux_yield;
    F.params = &p;

    scale = (1.0+redshift)*Wm2Hz_TO_mJy/(4.0*M_PI*pow(lumdist(redshift)*MPC_TO_METER,2.0));

    result = error = -1;
    bounds[0] = filters.get(filter_id).low();
    bounds[1] = filters.get(filter_id).high();
    
    gsl_integration_qags (&F, bounds[0], bounds[1], 0, SL_INT_PRECISION, SL_INT_SIZE, w, &result, &error); 

    result*=scale;
    fluxes[filter_id][params] = result;
  }
  
  return fluxes[filter_id][params];
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
