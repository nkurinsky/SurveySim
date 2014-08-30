#include "simulator.h"
#include "functions.h"
#include "constants.h"
#include "cosmo.h"

sprop::sprop(double z,double f[],double lum, double w){
  redshift = z;
  luminosity = lum;
  weight = w;

  bool do_colors(true);
  for (int i=0;i<3;i++){
    fluxes[i] = f[i];
    if(f[i] < 0.00000000E0)
      do_colors = false;
  }
  if(do_colors){
    c1 = get_color(f[2],f[0]);
    c2 = get_color(f[2],f[1]);
  }
  else{
    c1 = -99.0;
    c2 = -99.0;
  }
}

sprop::sprop(){
  redshift = -1;
  luminosity = -1;
  weight = 0;

  c1 = -99.0;
  c2 = -99.0;
}

products::products(int nz, int ns[]){
  dndz.resize(nz);
  dnds[0].resize(ns[0]);
  dnds[1].resize(ns[1]);
  dnds[2].resize(ns[2]);
}

clementsBins::clementsBins(){
  bnum[0] = 16;
  bnum[1] = 13;
  bnum[2] = 10;
  double btemp250[] = {34.0,38.0,42.0,47.0,52.0,56.0,63.0,73.0,83.0,93.0,125.0,175.0,250.0,350.0,450.0,800.0};
  double btemp350[] = {38.0,41.0,46.0,50.0,55.0,61.0,71.0,82.0,91.0,125.0,175.0,250.0,700.0};
  double btemp500[] = {48.0,52.0,56.0,60.0,63.0,70.0,82.0,95.0,150.0,250.0};

  double dbtemp250[] = {0.004,0.004,0.005,0.005,0.004,0.007,0.01,0.01,0.01,0.032,0.05,0.075,0.1,0.1,0.35,0.5};
  double dbtemp350[] = {0.003,0.005,0.004,0.005,0.006,0.01,0.011,0.009,0.034,0.05,0.075,0.45,0.5};
  double dbtemp500[] = {0.004,0.004,0.004,0.003,0.007,0.012,0.013,0.055,0.1,0.25};

  double Stemp250[] = {0.000213156,0.000281487,0.000361512,0.0004789,0.000616607,0.000742113,0.000996211,0.00143982,0.0019847,0.00263759,0.00552427,0.0128114,0.03125,0.072472,0.135841,0.572433};
  double Stemp350[] = {0.000281487,0.000340377,0.000453831,0.000559017,0.000709425,0.000919019,0.00134322,0.00192546,0.00249806,0.00552427,0.0128114,0.03125,0.409963};
  double Stemp500[] = {0.000504781,0.000616607,0.000742113,0.000881816,0.000996211,0.00129642,0.00192546,0.00278169,0.00871421,0.03125};
 
  b250.resize(bnum[0]);
  b350.resize(bnum[1]);
  b500.resize(bnum[2]);
  db250.resize(bnum[0]);
  db350.resize(bnum[1]);
  db500.resize(bnum[2]);
  S250.resize(bnum[0]);
  S350.resize(bnum[1]);
  S500.resize(bnum[2]);

  for (int i=0;i<bnum[0];i++){
    b250[i] = btemp250[i];
    db250[i] = dbtemp250[i];
    S250[i] = Stemp250[i];
  }
  
  for (int i=0;i<bnum[1];i++){
    b350[i] = btemp350[i];
    db350[i] = dbtemp350[i];
    S350[i] = Stemp350[i];
  }
  
  for (int i=0;i<bnum[2];i++){
    b500[i] = btemp500[i];
    db500[i] = dbtemp500[i];
    S500[i] = Stemp500[i];
  }
  
}

int simulator::binNum(int band, double flux){  
  switch(band){
  case 0:
    for (int i=0;i<dndsInfo.bnum[0];i++)
      if(flux < dndsInfo.b250[i])
	return i-1;
    break;
  case 1:
    for (int i=0;i<dndsInfo.bnum[1];i++)
      if(flux < dndsInfo.b350[i])
	return i-1;
    break;
  case 2:
    for (int i=0;i<dndsInfo.bnum[2];i++)      
      if(flux < dndsInfo.b500[i])
	return i-1;
    break;
  default:
    cout << "ERROR: Invalid band request, simulator::binNum\n";
  }
  
  return dndsInfo.bnum[band]-1;
}

simulator::simulator(string filterfile, string obsfile, string sedfile){

  seds.reset(NULL);

  color_exp = 0.0;
  printf("Initializing Observations:\n");
  observations.reset(new obs_lib(obsfile));
  printf("Initializing SED Library:\n");
  seds.reset(new sed_lib(sedfile));
  
  string filters[3];
  double f_lims[3];
  double errors[3];
  printf("Getting Filter Information\n");
  observations->info(filters,f_lims,errors);

  //initialize filters
  printf("Initializing Filters\n");
  if(seds->init_filter_lib(filterfile)){
    for(short i=0;i<3;i++){
      seds->load_filter(i,filters[i]);
      flux_limits[i]= f_lims[i];
      band_errs[i] = errors[i];
    }
  }

  //initialize observation portion of histograms
  //possibility here to integrate obs_lib into hist_lib
  printf("Initializing Diagnostic Class\n");
  diagnostic.reset(new hist_lib());
  double *c1,*c2;
  int osize = observations->get_snum();
  printf("\tGetting Colors\n");
  observations->get_all_colors(c1,c2);
  printf("\tInitializing Observations\n");
  diagnostic->init_obs(c1,c2,osize);
  last_output.chisqr=0;

  printf("Initializing Counts Containers\n");
  last_output.dnds[0].resize(dndsInfo.bnum[0]);
  last_output.dnds[1].resize(dndsInfo.bnum[1]);
  last_output.dnds[2].resize(dndsInfo.bnum[2]);
  printf("Simulator Initialized\n");
}

bool simulator::load_filter_lib(string file){
  return seds->init_filter_lib(file);
}

bool simulator::load_filter(short filt_id, string name, double error, double flim){
  if((filt_id < 3) and (filt_id >= 0)){
    band_errs[filt_id] = error;
    flux_limits[filt_id] = flim;
  }
  else{
    printf("ERROR: Invalid filt_id\n");
    return false;
  }
  return seds->load_filter(filt_id,name);
}

void simulator::set_size(double area,double dz,double zmin,int nz,int ns){
  this->area = (area > 0) ? ((area < 41254.0) ? area : 41253.0) : 1;
  this->dz = (dz > 0) ? dz : 0.1;
  this->zmin = (zmin > 0.001) ? zmin : 0.001;
  this->nz = (nz > 1) ? nz : 50;
  this->ns = (ns > 1) ? ns : 8; //currently unused

  last_output.dndz.resize(nz);
}

void simulator::set_color_exp(double val, double zcut){
  if (seds.get() != NULL)
    seds->set_color_evolution(val,zcut);
}

void simulator::set_lumfunct(lumfunct *lf){
  if(lf != NULL){
    this->lf = lf;
  }
  else
    cout << "ERROR: NULL Pointer Passed to Simulator" << endl;
}


void simulator::set_sed_lib(string sedfile){
  seds.reset(new sed_lib(sedfile));
}


void simulator::set_obs(string obsfile){
  observations.reset(new obs_lib(obsfile));
  double *c1,*c2;
  int osize = observations->get_snum();
  observations->get_all_colors(c1,c2);
  diagnostic->init_obs(c1,c2,osize);
}

products simulator::simulate(){
  sources.clear();
  products output(nz,dndsInfo.bnum);
  
  if(seds.get() != NULL){

    static int is,js,jsmin;
    static double tmpz,vol;

    int lnum = seds->get_lnum();
    double dl = seds->get_dl();
    static long nsrcs;
    double lums[lnum];
    double zarray[nz];    
    
    seds->get_lums(lums);
    
    //=========================================================================
    // for each L-z depending on the number of sources, sample the SED and get the appropriate fluxes
    //*************************************************************************
    
    //previous approach was round to nearest template value.
    // We interpolate here, in future will convolve filter
 
    static double flux_sim[3], flux_raw[3];
    static double noise[3];
    static sprop *temp_src;
    static int src_iter;
    static bool detected = true;
    static int dndsi;

    gsl_rng * r;
    const gsl_rng_type * T;    
    gsl_rng_default_seed = time(NULL);
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);    

    //NOTE templates are given in W/Hz
    for (is=0;is<nz;is++){
      zarray[is]=(is)*dz+zmin;

      tmpz=zarray[is]+dz/2.0;
      vol=(dvdz(tmpz,area)*dz);

      jsmin = 0;
      for (js=0;js<lnum;js++){
	flux_sim[0] = seds->get_filter_flux(lums[js],zarray[is],0);
	if(flux_sim[0]>=flux_limits[0]){
	  jsmin = (js > 0) ? (js-1) : js; //js-1 to allow for noise
	  js = lnum; //break out of loop
	};
      };
      //printf("Z: %f, Lmin: %f\n",zarray[is],lums[jsmin]);
      
      for (js=jsmin;js<lnum;js++){
	//source number is dn/(dldv)*Dv*Dl
	nsrcs = long(dl*vol*lf->get_nsrcs(zarray[is],lums[js]));
	for (int i=0;i<3;i++){
	  flux_raw[i] = seds->get_filter_flux(lums[js],zarray[is],i);
	}
	for (src_iter=0;src_iter<nsrcs;src_iter++){
	  detected = true;
	  for (int i=0;i<3;i++){
	    noise[i]=gsl_ran_gaussian(r,band_errs[i]);
	    flux_sim[i] = flux_raw[i]+noise[i];
	    if (flux_sim[i] < flux_limits[i]) //reject sources below flux limit
	      detected = false;
	    else{
	      dndsi = binNum(i,flux_sim[i]);
	      //printf("%i %i\n",i,dndsi);
	      if((dndsi >= 0) and (dndsi < dndsInfo.bnum[i]))
		output.dnds[i][dndsi]+=1.0;
	    }
	  }
	  
	  //check for detectability, if "Yes" add to list
	  if(detected){
	    temp_src = new sprop(zarray[is],flux_sim,lums[js],1.0);
	    sources.push_back(*temp_src);
	    output.dndz[is]+=1; 
	    delete temp_src;
	  }
	}
      }
    }
     
    //=========================================================================
    // generate diagnostic color-color plots
    //*************************************************************************
    
    int snum(sources.size());
    double * c1(new double[snum]);
    double * c2(new double[snum]);
    double * w (new double[snum]);
    for (int i=0;i<snum;i++){
      c1[i] = sources[i].c1;
      c2[i] = sources[i].c2;
      w[i] = sources[i].weight;
    }
    
    diagnostic->init_model(c1,c2,w,snum);
    output.chisqr=diagnostic->get_chisq();
    //note: this is not general, but temporarily implemented for publication
    output.dnds[0]*=dndsInfo.S250/dndsInfo.db250;
    output.dnds[1]*=dndsInfo.S350/dndsInfo.db350;
    output.dnds[2]*=dndsInfo.S500/dndsInfo.db500;
    output.dnds[0]/=area;
    output.dnds[1]/=area;
    output.dnds[2]/=area;
    last_output = output;

    //printf("%f %f %f\n",output.dnds[0][0],output.dnds[1][0],output.dnds[2][0]);
    
    delete[] c1;
    delete[] c2;
    delete[] w;
    
  }
  else{
    cout << "ERROR: NULL Model Library" << endl;
  }
  return output;
}


bool simulator::save(string outfile){
  try{
    printf("\tWriting Diagnostic\n");
    bool opened =  diagnostic->write_fits(outfile);
    printf("\tDiagnostic Written\n");
    
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
    printf("%s %lu\n","Sources Being Saved: ",size);
    
    valarray<double> f1(size),f2(size),f3(size),luminosity(size),redshift(size);
    for(unsigned long i=0;i<size;i++){
      f1[i] = sources[i].fluxes[0];
      f2[i] = sources[i].fluxes[1];
      f3[i] = sources[i].fluxes[2];
      redshift[i] = sources[i].redshift;
      luminosity[i] = sources[i].luminosity;
    }
    
    printf("1");
    
    static std::vector<string> colname(11,"");
    static std::vector<string> colunit(11,"-");
    static std::vector<string> colform(11,"f13.8");
    
    colname[0] = "F1";
    colname[1] = "F2";
    colname[2] = "F3";
    colname[3] = "Z";
    colname[4] = "Lum";
    colname[5] = "s1";
    colname[6] = "s2";
    colname[7] = "s3";
    colname[8] = "dNds1";
    colname[9] = "dNds2";
    colname[10] = "dNds3";
    
    colunit[0] = "Jy";
    colunit[1] = "Jy";
    colunit[2] = "Jy";
    colunit[4] = "W/Hz";
    colunit[5] = "mJy";
    colunit[6] = "mJy";
    colunit[7] = "mJy";
    colunit[8] = "Jy^1.5/sr";
    colunit[9] = "Jy^1.5/sr";
    colunit[10] = "Jy^1.5/sr";
    
    colform[4] = "e13.5";
    colform[8] = "e13.5";
    colform[9] = "e13.5";
    colform[10] = "e13.5";
    
    static string hname("Parameter Distributions");
    Table* newTable;
    try{
      newTable = pFits->addTable(hname,size,colname,colform,colunit,AsciiTbl);
    }
    catch(...){
      printf("Could not create table\n");
      exit(1);
    }
    
    printf("2");
    
    newTable->column(colname[0]).write(f1,1);
    newTable->column(colname[1]).write(f2,1);
    newTable->column(colname[2]).write(f3,1);
    newTable->column(colname[3]).write(redshift,1);
    newTable->column(colname[4]).write(luminosity,1);
    newTable->column(colname[5]).write(dndsInfo.b250,1);
    newTable->column(colname[6]).write(dndsInfo.b350,1);
    newTable->column(colname[7]).write(dndsInfo.b500,1);
    newTable->column(colname[8]).write(last_output.dnds[0],1);
    newTable->column(colname[9]).write(last_output.dnds[1],1);
    newTable->column(colname[10]).write(last_output.dnds[2],1);
    
    printf("3\n");

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
