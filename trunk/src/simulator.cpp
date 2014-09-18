#include "simulator.h"

simulator::simulator() : color_exp(0.0),
			 area(3.046174e-4), //default to 1sq degree
			 dz(0.1),
			 zmin(0.1), 
			 nz(59),
			 filterFile("/usr/local/surveysim/filters/filterlib.txt"), 
			 obsFile(""){
  
  char * ffile = getenv("FILTERFILE");
  if(ffile != NULL)
    filterFile = static_cast<string>(ffile);
  
  filters[0] = "F1";
  filters[1] = "F2";
  filters[2] = "F3";
  
  axes[0] = ColorF1F3;
  axes[1] = ColorF2F3;
}

simulator::simulator(string filterfile, string obsfile, string sedfile, axis_type axes[]) : simulator(){
 
  filterFile = filterfile;
  set_sed_lib(sedfile);
  this->axes[0] = axes[0];
  this->axes[1] = axes[1];
  set_obs(obsfile);
  
}

bool simulator::load_filter_lib(string file){
  filterFile = file;
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

void simulator::set_diagnostic_xaxis(axis_type option){
  axes[0] = option;
  reset_obs();
}

void simulator::set_diagnostic_yaxis(axis_type option){
  axes[1] = option;
  reset_obs();
}

void simulator::set_diagnostic_axes(axis_type xopt, axis_type yopt){
  axes[0] = xopt;
  axes[1] = yopt;
  reset_obs();
}

void simulator::set_size(double area,double dz,double zmin,int nz){
  this->area = (area > 0) ? ((area < 41254.0) ? area : 41253.0) : 1;
  this->dz = (dz > 0) ? dz : 0.1;
  this->zmin = (zmin > 0.001) ? zmin : 0.001;
  this->nz = (nz > 1) ? nz : 50;

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

void simulator::initialize_filters(){

  if((seds.get() != NULL) and (observations.get() != NULL)){
    
    if(seds->init_filter_lib(filterFile)){
      for(short i=0;i<3;i++)
	seds->load_filter(i,filters[i]);
    }
  }
  
}

void simulator::initialize_counts(){
  if(observations.get() != NULL){
    
    int osize = observations->get_snum();
    valarray<double> fluxes[3];

    for(int i=0;i<3;i++){
      fluxes[i].resize(osize,0.0);
      for(int j=0;j<osize;j++)
	fluxes[i][j] = observations->get_flux(j,i);
      counts[i].reset(new NumberCounts(fluxes[i],area,filters[i]));
      last_output.dnds[i].resize(counts[i]->bins().size());
    }

  }
}

void simulator::set_sed_lib(string sedfile){
  seds.reset(new sed_lib(sedfile));
  initialize_filters();
}

void simulator::set_obs(string obsfile){
  if(obsfile != ""){
    obsFile=obsfile;
    reset_obs();
    
    observations->info(filters,flux_limits,band_errs);
    initialize_filters();
    
    last_output.chisqr=0;  
    initialize_counts();
  }
  else
    printf("simulator::set_obs Error: tried to set obs_lib with empty file name\n");
}

void simulator::reset_obs(){
  if(obsFile != ""){
    diagnostic.reset(new hist_lib());
    observations.reset(new obs_lib(obsFile, axes));
    
    double *x,*y;
    int osize = observations->get_snum();
    observations->get_all_colors(x,y);
    diagnostic->init_obs(x,y,osize);
  }
}

products simulator::simulate(){
  sources.clear();
  int ns[] = {static_cast<int>(counts[0]->bins().size()),
	      static_cast<int>(counts[1]->bins().size()),
	      static_cast<int>(counts[2]->bins().size())};
  products output(nz,ns);
  
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
    static sprop *temp_src;
    static int src_iter;
    static bool detected = true;    

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
	    flux_sim[i] = rng.gaussian(flux_raw[i],band_errs[i],0.0,1e5);
	    if (flux_sim[i] < flux_limits[i]) //reject sources below flux limit
	      detected = false;
	  }
	  
	  //check for detectability, if "Yes" add to list
	  if(detected){
	    temp_src = new sprop(zarray[is],flux_sim,lums[js],1.0,axes);
	    sources.push_back(*temp_src);
	    output.dndz[is]++; 
	    delete temp_src;
	  }
	}
      }
    }
     
    //=========================================================================
    // generate diagnostic color-color plots
    //*************************************************************************
    
    int snum(sources.size());
    vector<double> c1(snum,0.0);
    vector<double> c2(snum,0.0);
    vector<double> w (snum,0.0);
    valarray<double> f1(0.0,snum);
    valarray<double> f2(0.0,snum);
    valarray<double> f3(0.0,snum);
    for (int i=0;i<snum;i++){
      c1[i] = sources[i].c1;
      c2[i] = sources[i].c2;
      w[i] = sources[i].weight;
      f1[i] = sources[i].fluxes[0];
      f2[i] = sources[i].fluxes[1];
      f3[i] = sources[i].fluxes[2];
    }
    
    diagnostic->init_model(c1.data(),c2.data(),w.data(),snum);
    output.chisqr=diagnostic->get_chisq();
    counts[0]->compute(f1,area,output.dnds[0]);
    counts[1]->compute(f2,area,output.dnds[1]);
    counts[2]->compute(f3,area,output.dnds[2]);
    last_output = output;
    
  }
  else{
    printf("ERROR: NULL Model Library\n");
  }
  return output;
}


bool simulator::save(string outfile){
  try{
    bool opened =  diagnostic->write_fits(outfile);
    
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
    
    static std::vector<string> colname(14,"");
    static std::vector<string> colunit(14,"-");
    static std::vector<string> colform(14,"f13.8");
    
    colname[0] = "F1";
    colname[1] = "F2";
    colname[2] = "F3";
    colname[3] = "Z";
    colname[4] = "Lum";
    colname[5] = "s1";
    colname[6] = "s2";
    colname[7] = "s3";
    colname[8] = "obs_dnds1";
    colname[9] = "obs_dnds2";
    colname[10] = "obs_dnds3";
    colname[11] = "mod_dnds1";
    colname[12] = "mod_dnds2";
    colname[13] = "mod_dnds3";

    
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
    colunit[11] = "Jy^1.5/sr";
    colunit[12] = "Jy^1.5/sr";
    colunit[13] = "Jy^1.5/sr";
    
    colform[4] = "e13.5";
    colform[8] = "e13.5";
    colform[9] = "e13.5";
    colform[10] = "e13.5";
    colform[11] = "e13.5";
    colform[12] = "e13.5";
    colform[13] = "e13.5";
    
    static string hname("Parameter Distributions");
    Table* newTable;
    try{
      newTable = pFits->addTable(hname,size,colname,colform,colunit,AsciiTbl);
    }
    catch(...){
      printf("Could not create table\n");
      exit(1);
    }
    
    newTable->column(colname[0]).write(f1,1);
    newTable->column(colname[1]).write(f2,1);
    newTable->column(colname[2]).write(f3,1);
    newTable->column(colname[3]).write(redshift,1);
    newTable->column(colname[4]).write(luminosity,1);

    valarray<double> counts1(pow(10,counts[0]->bins()));
    valarray<double> counts2(pow(10,counts[1]->bins()));
    valarray<double> counts3(pow(10,counts[2]->bins()));

    newTable->column(colname[5]).write(counts1,1);
    newTable->column(colname[6]).write(counts2,1);
    newTable->column(colname[7]).write(counts3,1);

    newTable->column(colname[8]).write(counts[0]->counts(),1);
    newTable->column(colname[9]).write(counts[1]->counts(),1);
    newTable->column(colname[10]).write(counts[2]->counts(),1);

    newTable->column(colname[11]).write(last_output.dnds[0],1);
    newTable->column(colname[12]).write(last_output.dnds[1],1);
    newTable->column(colname[13]).write(last_output.dnds[2],1);

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
