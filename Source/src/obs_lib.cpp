#include "obs_lib.h"

obs::obs(){
  c1 = -99;
  c2 = -99;
}

obs::obs(double * f,double *bands){
  for (int i=0;i<3;i++){
    fluxes[i] = f[i];
  }
  
  c1 = get_color(fluxes[2],fluxes[0]);
  c2 = get_color(fluxes[2],fluxes[1]);
}

double obs::get_flux(int band){
  switch(band){
  case 0:
  case 1:
  case 2:
    return fluxes[band];
    break;
  default:
    return -1;
  }
}

void obs::get_colors(double &c1,double &c2){
  c1 = this->c1;
  c2 = this->c2;
}

obs::~obs(){
  if(interp_init){
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
}

obs_lib::obs_lib(string fitsfile){
  //initialize FITS input
  std::auto_ptr<FITS> pInfile(new FITS(fitsfile,Read));
  
  //reading primary header from fits file
  HDU& header = pInfile->pHDU();
  int hdunum(1);
  
  try{header.readKey("FHDU",hdunum);}
  catch(HDU::NoSuchKeyword){
    printf("FHDU not set, defaulting to HDU %i\n",extnum);}
  
  ExtHDU& table;
  try{table = pInfile->extension(extnum);}
  catch(FITS::NoSuchHDU){
    printf("HDU %i does not exist inf %s\n",extnum,fitsfile.c_str());
    exit(1);
  }

  const string ti[3] = {"1","2","3"};
  string column[6];
  for (int i=0;i<3;i++){
    try{table.readKey("F"+ti[i]+"COL",column[i]);}
    catch(HDU::NoSuchKeyword){
      printf("Keyword F%iCOL missing from header (%s)\n",i,fitsfile.c_str());
      exit(1);
    }
    try{table.readKey("EF"+ti[i]+"COL",column[i+3]);}
    catch(HDU::NoSuchKeyword){
      printf("Keyword EF%iCOL missing from header (%s)\n",i,fitsfile.c_str());
      exit(1);
    }
    try{table.readKey("F"+ti[i]+"MIN",flim[i]);}
    catch(HDU::NoSuchKeyword){
      printf("Keyword F%iLIM missing from header (%s)\n",i,fitsfile.c_str());
      exit(1);
    }
    try{table.readKey("F"+ti[i]+"FILT",filter[i]);}
    catch(HDU::NoSuchKeyword){
      printf("Keyword F%iFILT missing from header (%s)\n",i,fitsfile.c_str());
      exit(1);
    }
  }
  
  unsigned long tablesize(table.rows());
  valarray<double> col(tablesize)[6];
  string unit;
  for(int i=0;i<6;i++){
    try{
      unit = table.column(column[i]).unit();
      table.column(column[i]).read(col[i],1,tablesize);
      if(unit == "Jy")
	col[i] *= 1e3; //convert to mJy
    }
    catch(Table::NoSuchColum){
      printf("Column %s does not exist in %s\n",column[i],fitsfile.c_str());
      exit(1);
    }
  }

  //errors are mean of observation errors
  for(int i=0;i<3;i++){
    ferr[i] = col[i+3].sum()/tablesize;
  }

  observations.resize(tablesize);
  double fluxes[3];
  
  for (unsigned int i=0;i<observations.size();i++){
    for (int j=0;j<3;j++)
      fluxes[j]=col[j][i];
    observations[i].reset(new obs(fluxes));
  }
  
}
  
double obs_lib::get_flux(int i,int band){
  if(i < int(observations.size()))
    return observations[i]->get_flux(band);
  else
    return -1;
}

void obs_lib::get_colors(int i,double &c1,double &c2){
  if(i < int(observations.size()))
    observations[i]->get_colors(c1,c2);
  else{
    c1 = -1;
    c2 = -1;
  }
}

void obs_lib::get_all_colors(double* &c1,double* &c2){
  c1 = new double[observations.size()];
  c2 = new double[observations.size()];
  for (unsigned int i=0;i<observations.size();i++){
    observations[i]->get_colors(c1[i],c2[i]);
  }
}

void obs_lib::info(string filters[],double flims[],double ferrs[]){
  if(filters == NULL)
    filters = new double[3];
  if(flims == NULL)
    flims = new double[3];
  if(ferrs == NULL)
    ferrs = new double[3];
  
  for(int i=0;i<3;i++){
    filters[i] = this->filter[i];
    flims[i] = this->flim[i];
    ferrs[i] = this->ferr[i];
  }
}

obs_lib::~obs_lib(){

}
