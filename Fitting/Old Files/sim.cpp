#include "sim.h"

sprop::sprop(double m,double z,double b[],double f[],double lum){
  mod_num = m;
  redshift = z;
  luminosity = lum;
  bool do_colors(true);
  for (int i=0;i<3;i++){
    bands[i] = b[i];
    fluxes[i] = f[i];
    if((b[i] < 0.0000000E0) or (f[i] < 0.00000000E0))
      do_colors = false;
  }
  if(do_colors){
    c1 = log10(f[2]/f[0])/log10(b[2]/b[0]);
    c2 = log10(f[2]/f[1])/log10(b[2]/b[1]);
  }
  else{
    c1 = -99.0;
    c2 = -99.0;
  }
}

sprop::sprop(){
  mod_num = -1;
  redshift = -1;
  luminosity = -1;
  for (int i=0;i<3;i++){
    fluxes[i] = -1;
    bands[i] = -1;
  }
  c1 = -99.0;
  c2 = -99.0;
}

sim::sim(model_lib* mods){
  lib = mods;
  simulated = false;
  initialized = false;
  pnum = lib->get_pnum();
  init_rand();
  this->mods = new distribution*[pnum];
}

void sim::initialize_dists(double mean[],double sigma[],int size,double min[],double max[]){
  static double p_range[2],z_range[2];
  z_range[0] = min[0];
  z_range[1] = max[0];

  if(initialized)
    delete z;
  if(mean[0] <= 0){
    z = new distribution(r,size,z_range);}
  else if (sigma[0] > 0){
    z = new distribution(r,size,mean[0],sigma[0],z_range);}
  else{
    z = new distribution(size,mean[0]);}
  
  for (int i=1;i <= pnum;i++){
    if(initialized)
      delete mods[i-1];
    p_range[0] = min[i];
    p_range[1] = max[i];
    if(mean[i] <= 0){
      mods[i-1] = new distribution(r,size,p_range);}
    else if (sigma[i] > 0){
      mods[i-1] = new distribution(r,size,mean[i],sigma[i],p_range);}
    else{
      mods[i-1] = new distribution(size,mean[i]);}
  }

  initialized = true;
}

void sim::init_rand(){
  gsl_rng_default_seed = time(NULL);
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
}

void sim::add_dists(double mean[],double sigma[],int size,double min[],double max[]){
  static double p_range[2],z_range[2];
  z_range[0] = min[0];
  z_range[1] = max[0];

  if(mean[0] <= 0){
    z->add(r,size,z_range);}
  else if (sigma[0] > 0){
    z->add(r,size,mean[0],sigma[0],z_range);}
  else{
    z->add(size,mean[0]);}
  
  for (int i=1;i <= pnum;i++){
    p_range[0] = min[i];
    p_range[1] = max[i];
    if(mean[i] <= 0){
      mods[i-1]->add(r,size,p_range);}
    else if (sigma[i] > 0){
      mods[i-1]->add(r,size,mean[i],sigma[i],p_range);}
    else{
      mods[i-1]->add(size,mean[i]);}
  }
}

sprop sim::get(int i){
  if (simulated){
    return sources[i];
  }
  else{
    static sprop retval;
    return retval;
  }
}

void sim::get(int i,double &rz,double m[]){
  rz = z->get(i);
  for (int j=0;j<pnum;j++){
    m[j] = mods[j]->get(i);
  }
}

void sim::get(int i,double &rz,double &m){
  int tmods[pnum];
  rz = z->get(i);
  for (int j=0;j<pnum;j++){
    tmods[j] = mods[j]->get(i);
  }
  m = lib->index(tmods);
}

int sim::size(){
  for (int i=0;i<pnum;i++){
    if(z->size() != mods[i]->size()) 
      return -1;
  }
  return z->size();
}

void sim::get_color(int i,double &c1,double &c2){
  if(simulated){
    if(i < int(sources.size())){
      c1 = sources[i].c1;
      c2 = sources[i].c2;
      return;
    }
    else{
      cout << "i out of range" << endl;
    }
  }
  
  c1 = -99.0;
  c2 = -99.0;
  return;
}

void sim::get_all_colors(double *&c1,double *&c2){
  if(simulated){
    int snum(sources.size());
    c1 = new double[snum];
    c2 = new double[snum];
    for (int i=0;i<snum;i++){
      c1[i] = sources[i].c1;
      c2[i] = sources[i].c2;
    }
  }
  else{
    c1 = NULL;
    c2 = NULL;
  }
}

bool sim::simulate(double bands[]){
  if (size() == -1)
    return false;
  
  static sprop *temp;
  static double rz;
  double mod_nums[pnum];
  static double fluxes[3];

  for (int i=0;i<size();i++){
    get(i,rz,mod_nums);
    for (int j=0;j<3;j++){
      fluxes[j] = lib->get_flux(mod_nums,bands[j],rz);
    }
    temp = new sprop(lib->index(mod_nums),rz,bands,fluxes,-1);
    sources.push_back(*temp);
    delete temp;
  }
  simulated = true;
  return simulated;
}

bool sim::simulate(double bands[],double bfluxerr[],double flim[],lumfunct &lums){
  if (size() == -1)
    return false;
  
  sources.clear();

  static sprop *temp;
  static double rz;
  double mod_nums[pnum];
  static double fluxes[3];
  static double f350,scale;
  static bool repeat(true);

  static double nrange[2];
  double ** noise = new double*[3];
  for (int i=0;i<3;i++){
    nrange[0] = bfluxerr[i]*-5.0;
    nrange[1] = bfluxerr[i]*5.0;
    noise[i] = gauss_random(r,nrange,0.0,bfluxerr[i],size());
  }

  static double ftemp[3],ltemp;
  static int itnum,errnum;

  errnum = 0;
  for (int i=0;i<size();i++){
    get(i,rz,mod_nums);
    itnum = 0;
    do{
      f350 = lums.get_amp(rz,ltemp);
      scale = f350/lib->get_flux(mod_nums,350.0e-6,rz);
      for (int j=0;j<3;j++){
	fluxes[j] = lib->get_flux(mod_nums,bands[j],rz);
	fluxes[j] *= scale;
	ftemp[j] = fluxes[j];
	fluxes[j] += noise[j][i];
      }
      itnum++;
      if(itnum == 10){
	errnum++;
	repeat = false;
	fluxes[2] = -1.0;
      }
      else if((ftemp[0] < flim[0]) or (ftemp[1] < flim[1]) or (ftemp[2] < flim[2]))
	repeat = true;
      else
	repeat = false;
    }while (repeat);
    temp = new sprop(lib->index(mod_nums),rz,bands,fluxes,ltemp);
    sources.push_back(*temp);
    delete temp;
  }

  for (int i=0;i<3;i++)
    delete[] noise[i];
  delete[] noise;

  //cout << "Number of Undetectable Sources: " << errnum << endl;

  simulated = true;
  return simulated;
}

sim::~sim(){
  delete z;
  for (int i=0;i<pnum;i++)
    delete mods[i];
  delete[] mods;
  gsl_rng_free(r);
}

ostream &operator<<(ostream &out,sprop c){
  out << "Model:" << '\t' << c.mod_num << '\t';
  out << "Redshift:" << '\t' << c.redshift << endl;
  out << "Bands:" << '\t';
  for (int i=0;i<3;i++)
    out << c.bands[i] << '\t';
  out << endl;
  out << "Fluxes:" << '\t';
  for (int i=0;i<3;i++)
    out << c.fluxes[i] << '\t';
  out << endl;
  out << "Colors:" << '\t' << c.c1 << '\t' << c.c2 << endl;
  return out;
}
