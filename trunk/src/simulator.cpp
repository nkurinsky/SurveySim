#include "simulator.h"

//Lets assume we always have configuration so no need to have a "default" mode
simulator::simulator(const Configuration &config){
  filters[0] = "F1";
  filters[1] = "F2";
  filters[2] = "F3";
  obsFile="";  
  color_exp=0.0; //default color evolution

  configure(config);
}

void simulator::configure(const Configuration &config){
  logflag=config.oprint;

  modelFile = config.modfile;
  double tarea(config.areaSteradian());

  area = (tarea > 0) ? ((tarea <= 12.56637) ? tarea : 12.56637) : 3.046174e-4;
  dz = (config.dz > 0) ? config.dz : 0.1;
  zmin = (config.zmin > 0.001) ? config.zmin : 0.001;
  nz = (config.nz > 1) ? config.nz : 50;
  last_output.dndz.resize(nz);
  last_output.chisqr=0;    
    
  seds.reset(new sed_lib(config.sedfile, nz, zmin, dz));
  seds->load_filters(modelFile,logflag);
  seds->get_filter_info(filters,flux_limits,band_errs);
  axes[0] = config.axes[0];
  axes[1] = config.axes[1];
  fagns.reset(new agn_frac(seds->get_tnum()));
  
  simflag=config.simflag;

  if( not simflag ){
    if(config.obsfile != ""){
      obsFile=config.obsfile;
      diagnostic.reset(new hist_lib());

      observations.reset(new obs_lib(obsFile, axes, flux_limits));
      
      vector<double> x,y;
      int osize = observations->get_snum();
      observations->get_all_colors(x,y);
      diagnostic->init_obs(x.data(),y.data(),osize);
    
      valarray<double> fluxes[3];
      for(int i=0;i<3;i++){
	fluxes[i].resize(osize,0.0);
	for(int j=0;j<osize;j++)
	  fluxes[i][j] = observations->get_flux(j,i);
	counts[i].reset(new NumberCounts(fluxes[i],area,filters[i]));
	last_output.dnds[i].resize(counts[i]->bins().size());
      }
    }
    else{
      LOG_CRITICAL(printf("simulator::configure Error: tried to set obs_lib with empty file name\n"));
    }
  }
  else{ //no observations, will be initialized at first simulation
    for(int i=0;i<3;i++){
      counts[i].reset();
      last_output.dnds[i].resize(0);
    }
  }
}

void simulator::set_color_exp(double val, double zcut){
  if (seds)
    seds->set_color_evolution(val,zcut);
}

void simulator::set_lumfunct(lumfunct *lf){
  if(lf != NULL){
    this->lf = lf;
  }
  else
    cout << "ERROR: NULL Pointer Passed to Simulator" << endl;
}

long simulator::num_sources(double z, double l, double dl){
  
  //minimum and maximum
  double zval[2]={z,z+dz};
  double lval[2]={l-dl/2.0,l+dl/2.0};

  //source number is integral(dn/(dldv)*(dv/dz),dl,dz)
  double nsrcs = 0.0;
  for(int i=0;i<2;i++){
    for(int j=0;j<2;j++){
      nsrcs += lf->get_phi(zval[j],lval[i])*dvdz(zval[j],area);
    }
  }
  //Multiply by dlogl and dz
  nsrcs*=0.25*dl*dz;
  
  long retval=static_cast<long>(floor(nsrcs));
  double p_extra_source = nsrcs-floor(nsrcs);
  
  if ((p_extra_source > 0.0) and (p_extra_source < 1.0) and (p_extra_source > rng.flat(0,1))){
    retval++;
  }
  else if (p_extra_source < 0.0){
    LOG_CRITICAL(printf("ERROR: val-floor less than 0!\n"));
    exit(1);
  }
  else if (p_extra_source > 1.0){
    LOG_CRITICAL(printf("ERROR: val-floor greater than 1!\n"));
    exit(1);
  }
  
  return retval;
}

void simulator::initial_simulation(){
  
  int is,js,jsmin;
  int lnum = seds->get_lnum();
  double dl = seds->get_dl();
  long nsrcs;
  double lums[lnum];
  double zarray[nz];
  short sedtype; // this holds the id of the SED type 
  
  seds->get_lums(lums);

  unique_ptr<int[]> jsmin_type(new int[seds->get_tnum()]);
  
  //=========================================================================
  // for each L-z depending on the number of sources, sample the SED and get 
  // the appropriate fluxes
  //*************************************************************************
  
  static double flux_sim[3], flux_raw[3];
  static sprop *temp_src;
  static int src_iter;
  static bool detected = true;    

  static double tL,tZ;
  
  //NOTE templates are given in W/Hz
  for (is=0;is<nz;is++){

    zarray[is]=(is)*dz+zmin;    

    jsmin=lnum;
    // here we determine minimum luminosity to attempt to simulate
    for(int stype=0;stype<seds->get_tnum();stype++){
      jsmin_type[stype] = 0;
      for (js=0;js<lnum;js++){
	flux_sim[0] = seds->get_filter_flux(lums[js],zarray[is],stype,0);
	if(flux_sim[0]>=flux_limits[0]){
	  jsmin_type[stype] = (js > 0) ? (js-1) : js; //js-1 to allow for noise
	  js = lnum; //break out of loop
	}
      }
      if(jsmin_type[stype] < jsmin)
	jsmin=jsmin_type[stype];
    }

    for (js=jsmin;js<lnum;js++){
      nsrcs = num_sources(zarray[is],lums[js],dl);
      
      for (src_iter=0;src_iter<nsrcs;src_iter++){
	detected = true;
	//determine sedtype (agn type)
	sedtype=fagns->get_sedtype(lums[js],zarray[is]);
	
	//get uniformly distributed luminosity and redshift
	if(js == 0)
	  tL=rng.flat(lums[js],lums[js]+dl/2.0);
	else if (js == lnum-1)
	  tL=rng.flat(lums[js]-dl/2.0,lums[js]);
	else
	  tL=rng.flat(lums[js]-dl/2.0,lums[js]+dl/2.0);
	
	tZ=rng.flat(zarray[is],zarray[is]+dz);

	for (int i=0;i<3;i++){
	  flux_raw[i] = seds->get_filter_flux(tL,tZ,sedtype,i);
	}

	for (int i=0;i<3;i++){
	  flux_sim[i] = rng.gaussian(flux_raw[i],band_errs[i],0.0,1e5);
	  if (flux_sim[i] < flux_limits[i]) //reject sources below flux limit
	    detected = false;
	}
	
	//check for detectability, if "Yes" add to list
	if(detected){
	  temp_src = new sprop(tZ,flux_sim,tL,axes,sedtype);
	  sources.push_back(*temp_src);
	  delete temp_src;
	}
      }
    }
  }
  
  //=========================================================================
  // generate diagnostic plots
  //*************************************************************************
  
  int snum(sources.size());
  valarray<double> f1(0.0,snum);
  valarray<double> f2(0.0,snum);
  valarray<double> f3(0.0,snum);
  for (int i=0;i<snum;i++){
    f1[i] = sources[i].fluxes[0];
    f2[i] = sources[i].fluxes[1];
    f3[i] = sources[i].fluxes[2];
  }

  for(int i=0;i<3;i++){ //for each band
    if(not counts[i]){ //if counts not initialized
      //initialize counts based on simulated fluxes
      switch(i){
      case 0:
	counts[i].reset(new NumberCounts(f1,area,filters[i]));
	break;
      case 1:
	counts[i].reset(new NumberCounts(f2,area,filters[i]));
	break;
      case 2:
	counts[i].reset(new NumberCounts(f3,area,filters[i]));
	break;
      }
      //resize output
      last_output.dnds[i].resize(counts[i]->bins().size());
    }
  }
}

products simulator::simulate(){
  //the "sources" structure holds the simulated sources
  sources.clear();
  
  if(seds.get() == NULL){
    LOG_CRITICAL(printf("ERROR: NULL Model Library\n"));
    return last_output;
  }

  if((not counts[0]) or (not counts[1]) or (not counts[2])){
    LOG_DEBUG(printf("Initializing Counts Structure\n"));
    initial_simulation();
  }

  int ns[] = {static_cast<int>(counts[0]->bins().size()),
	      static_cast<int>(counts[1]->bins().size()),
	      static_cast<int>(counts[2]->bins().size())};
  products output(nz,ns);
  
  static int is,js,jsmin;
  
  int lnum = seds->get_lnum();
  double dl = seds->get_dl();
  static long nsrcs;
  double lums[lnum];
  double zarray[nz];
  short sedtype; // this holds the id of the SED type 
  
  seds->get_lums(lums);

  unique_ptr<int[]> jsmin_type(new int[seds->get_tnum()]);
  
  //=========================================================================
  // for each L-z depending on the number of sources, sample the SED and get 
  // the appropriate fluxes
  //*************************************************************************
  
  static double flux_sim[3], flux_raw[3];
  static sprop *temp_src;
  static int src_iter;
  static bool detected = true;    

  static double tL,tZ;

  unsigned long srcTotal=0;
  unsigned long detTotal=0;
  //NOTE templates are given in W/Hz
  for (is=0;is<nz;is++){

    zarray[is]=(is)*dz+zmin;    

    jsmin=lnum;
    // here we determine minimum luminosity to attempt to simulate
    for(int stype=0;stype<seds->get_tnum();stype++){
      jsmin_type[stype] = 0;
      for (js=0;js<lnum;js++){
	flux_sim[0] = seds->get_filter_flux(lums[js],zarray[is],stype,0);
	if(flux_sim[0]>=flux_limits[0]){
	  jsmin_type[stype] = (js > 0) ? (js-1) : js; //js-1 to allow for noise
	  js = lnum; //break out of loop
	}
      }
      if(jsmin_type[stype] < jsmin)
	jsmin=jsmin_type[stype];
    }
    
    for (js=jsmin;js<lnum;js++){
      nsrcs = num_sources(zarray[is],lums[js],dl);
      
      for (src_iter=0;src_iter<nsrcs;src_iter++){
	srcTotal++;
	if(((srcTotal % 10000) == 0) and (srcTotal != 0)){
	  //printf("%lu\t%lu\n",srcTotal,detTotal);
	  //trying to avoid infinite loops
	  if(detTotal > 10*observations->get_snum()){
	    printf("Detected sources %lu (%i times obsnum) too high, aborting\n",detTotal,10);
	    output.chisqr=100;
	    return output;
	  }
	  else if (srcTotal > 1e6){
	    printf("Simulated sources %lu too high (>1e6), aborting\n",srcTotal);
            output.chisqr=100;
            return output; 
	  }
	}
	
	detected = true;
	//determine sedtype (agn type)
	sedtype=fagns->get_sedtype(lums[js],zarray[is]);
	
	//get uniformly distributed luminosity and redshift
	if(js == 0)
	  tL=rng.flat(lums[js],lums[js]+dl/2.0);
	else if (js == lnum-1)
	  tL=rng.flat(lums[js]-dl/2.0,lums[js]);
	else
	  tL=rng.flat(lums[js]-dl/2.0,lums[js]+dl/2.0);
	
	tZ=rng.flat(zarray[is],zarray[is]+dz);

	for (int i=0;i<3;i++){
	  flux_raw[i] = seds->get_filter_flux(tL,tZ,sedtype,i);
	}

	for (int i=0;i<3;i++){
	  flux_sim[i] = rng.gaussian(flux_raw[i],band_errs[i],0.0,1e5);
	  if (flux_sim[i] < flux_limits[i]) //reject sources below flux limit
	    detected = false;
	}
	
	//check for detectability, if "Yes" add to list
	if(detected){
	  detTotal++;
	  temp_src = new sprop(tZ,flux_sim,tL,axes,sedtype);
	  sources.push_back(*temp_src);
	  output.dndz[is]++; 
	  delete temp_src;
	}
      }
    }
  }
  
  //=========================================================================
  // generate diagnostic plots
  //*************************************************************************
  
  int snum(sources.size());
  
  valarray<double> f1(0.0,snum);
  valarray<double> f2(0.0,snum);
  valarray<double> f3(0.0,snum);
  for (int i=0;i<snum;i++){
    f1[i] = sources[i].fluxes[0];
    f2[i] = sources[i].fluxes[1];
    f3[i] = sources[i].fluxes[2];
  }

  if (not simflag){ //comparing to observations; get diagnostics
    vector<double> c1(snum,0.0);
    vector<double> c2(snum,0.0);
    for (int i=0;i<snum;i++){
      c1[i] = sources[i].c1;
      c2[i] = sources[i].c2;
    }
    diagnostic->init_model(c1.data(),c2.data(),snum);
    output.chisqr=diagnostic->get_chisq();
  }

  //compute counts
  counts[0]->compute(f1,area,output.dnds[0]);
  counts[1]->compute(f2,area,output.dnds[1]);
  counts[2]->compute(f3,area,output.dnds[2]);

  //store last computed modeled counts
  last_output = output;
  
  return output;
}


bool simulator::save(string outfile){

  try{
    bool opened = true;
    if(not simflag){ //if there are observations, save diagnostics
      opened = diagnostic->write_fits(outfile);
    }
    else
      outfile="!"+outfile;
    
    if(not opened)
      return opened;
    
    using namespace CCfits;
    //FITS::setVerboseMode(true);
    std::auto_ptr<FITS> pFits(0);
    
    try{
      pFits.reset(new FITS(outfile,Write));
    }
    catch (CCfits::FITS::CantOpen){
      cerr << "Can't open " << outfile << endl;
      opened = false;       
    }
    
    double lpars[LUMPARS];
    lf->get_params(lpars);
    
    pFits->pHDU().addKey("PHI0",lpars[0],"Normalization"); 
    pFits->pHDU().addKey("L0",lpars[1],"Knee location z=0"); 
    pFits->pHDU().addKey("ALPHA",lpars[2],"upper slope"); 
    pFits->pHDU().addKey("BETA",lpars[3],"lower slope"); 
    pFits->pHDU().addKey("P",lpars[4],"Norm evolution"); 
    pFits->pHDU().addKey("Q",lpars[5],"Knee evolution"); 
    pFits->pHDU().addKey("ZCUT",lpars[6],"Z Evolution Cutoff");
    pFits->pHDU().addKey("CEXP",color_exp,"Color Evolution Param");
    
    unsigned long size = sources.size();
    LOG_INFO(printf("%s %lu\n","Sources Being Saved: ",size));
    
    valarray<double> f1(size),f2(size),f3(size),luminosity(size),redshift(size),sedtype(size);
    for(unsigned long i=0;i<size;i++){
      f1[i] = sources[i].fluxes[0];
      f2[i] = sources[i].fluxes[1];
      f3[i] = sources[i].fluxes[2];
      redshift[i] = sources[i].redshift;
      luminosity[i] = sources[i].luminosity;
      sedtype[i] = sources[i].sedtype;
    }
    
    //if simulation-only mode we only need 11columns, otherwise add observed counts
    int colnum=15;
    if(simflag)
      colnum-=3;

    static std::vector<string> colname(colnum,"");
    static std::vector<string> colunit(colnum,"-");
    static std::vector<string> colform(colnum,"e13.5");

    colname[0] = "F1";
    colname[1] = "F2";
    colname[2] = "F3";
    colname[3] = "Z";
    colname[4] = "Lum";
    colname[5] = "Type";
    colname[6] = "s1";
    colname[7] = "s2";
    colname[8] = "s3";
    colname[9] = "mod_dnds1";
    colname[10] = "mod_dnds2";
    colname[11] = "mod_dnds3";
    
    colunit[0] = "Jy";
    colunit[1] = "Jy";
    colunit[2] = "Jy";
    colunit[4] = "W/Hz";
    colunit[6] = "mJy";
    colunit[7] = "mJy";
    colunit[8] = "mJy";
    colunit[9] = "Jy^1.5/sr";
    colunit[10] = "Jy^1.5/sr";
    colunit[11] = "Jy^1.5/sr";

    if(!simflag) {
      colname[12] = "obs_dnds1";
      colname[13] = "obs_dnds2";
      colname[14] = "obs_dnds3";

      colunit[12] = "Jy^1.5/sr";
      colunit[13] = "Jy^1.5/sr";
      colunit[14] = "Jy^1.5/sr";
    }
    
    
    static string hname("Parameter Distributions");
    Table* newTable;
    try{
      newTable = pFits->addTable(hname,size,colname,colform,colunit,AsciiTbl);
    }
    catch(...){
      printf("Could not create table\n");
      exit(1);
    }
    
    try{
      newTable->column(colname[0]).write(f1,1);
      newTable->column(colname[1]).write(f2,1);
      newTable->column(colname[2]).write(f3,1);
      newTable->column(colname[3]).write(redshift,1);
      newTable->column(colname[4]).write(luminosity,1);
      newTable->column(colname[5]).write(sedtype,1);

      valarray<double> counts1(pow(10,counts[0]->bins()));
      valarray<double> counts2(pow(10,counts[1]->bins()));
      valarray<double> counts3(pow(10,counts[2]->bins()));

      newTable->column(colname[6]).write(counts1,1);
      newTable->column(colname[7]).write(counts2,1);
      newTable->column(colname[8]).write(counts3,1);
      
      newTable->column(colname[9]).write(last_output.dnds[0],1);
      newTable->column(colname[10]).write(last_output.dnds[1],1);
      newTable->column(colname[11]).write(last_output.dnds[2],1);

      if(not simflag){
	newTable->column(colname[12]).write(counts[0]->counts(),1);
	newTable->column(colname[13]).write(counts[1]->counts(),1);
	newTable->column(colname[14]).write(counts[2]->counts(),1);
      }
    }
    catch(FitsException &except){
      printf("Caught Save Error: Column Write -- ");
      printf("%s\n",except.message().c_str());
      exit(1);
    }
    catch(...){
      printf("Caught Save Error: Column Write\n");
      exit(1);
    }

  }
  catch(...){
    printf("Caught Save Error\n");
    exit(1);
  }
  
  return true;
}
  
simulator::~simulator(){

}

void simulator::reset(){
  sources.clear();
}
