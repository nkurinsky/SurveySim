#include "simulator_utils.h"

sprop::sprop(){
  redshift = -1;
  luminosity = -1;
  weight = 0;
  
  c1 = -99.0;
  c2 = -99.0;
}

sprop::sprop(double z,double f[],double lum, double w, axis_type opts[]) : sprop(){
  redshift = z;
  luminosity = lum;
  weight = w;

  bool do_metrics(true);
  for (int i=0;i<3;i++){
    fluxes[i] = f[i];
    if(f[i] < 0.00000000E0)
      do_metrics = false;
  }
  if(do_metrics){
    c1 = metric_value(f[0],f[1],f[2],opts[0]);
    c2 = metric_value(f[0],f[1],f[2],opts[1]);
  }
  
}


products::products(int nz, int ns[]) : chisqr(-1){
  dndz.resize(nz);
  dnds[0].resize(ns[0]);
  dnds[1].resize(ns[1]);
  dnds[2].resize(ns[2]);
}
