#include "obs_lib.h"

obs::obs(){
  c1 = -99;
  c2 = -99;
}

obs::obs(double * f, axis_type axes[]){
  for (int i=0;i<3;i++){
    fluxes[i] = f[i];
  }
  
  c1 = metric_value(fluxes[0],fluxes[1],fluxes[2],axes[0]);
  c2 = metric_value(fluxes[0],fluxes[1],fluxes[2],axes[1]);
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

}

obs_lib::obs_lib(string fitsfile, axis_type axes[], double flim[]){
  //initialize FITS input
  std::unique_ptr<FITS> pInfile;
  
  try{
    pInfile.reset(new FITS(fitsfile,Read));
  }
  catch(...){
    printf("Could not open file \"%s\"\n",fitsfile.c_str());
    exit(1);
  }

  //reading primary header from fits file
  HDU& header = pInfile->pHDU();
  int hdunum(1);
  
  try{header.readKey("FHDU",hdunum);}
  catch(HDU::NoSuchKeyword){
    printf("FHDU not set, defaulting to HDU %i\n",hdunum);}

  try{
  
    ExtHDU& table = pInfile->extension(hdunum);

    char cnum[2];
    string num;
    string column[6];
    for (int i=0;i<3;i++){
      sprintf(cnum,"%i",i+1);
      num = cnum;
      try{table.readKey("F"+num+"COL",column[i]);}
      catch(HDU::NoSuchKeyword){
	printf("Keyword F%sCOL missing from header (%s)\n",num.c_str(),fitsfile.c_str());
	exit(1);
      }
      try{table.readKey("EF"+num+"COL",column[i+3]);}
      catch(HDU::NoSuchKeyword){
	printf("Keyword EF%sCOL missing from header (%s)\n",num.c_str(),fitsfile.c_str());
	exit(1);
      }
    }
    
    unsigned long tablesize(table.rows());
    valarray<valarray<double> > col(valarray<double>(tablesize),6);
    string unit;
    zp[0]=8.9;
    zp[1]=8.9;
    zp[2]=8.9;

    for(int i=0;i<6;i++){
      try{
	unit = table.column(column[i]).unit();
	table.column(column[i]).read(col[i],1,tablesize);
	//	std::cout << unit << std::endl;
	if(unit == "Jy")
	  col[i] *= 1e3; //convert to mJy
	if(unit == "uJy")
	  col[i] /= 1e3; //convert to mJy
	if(unit == "mag"){ //works for AB magnitudes only, should consider including Vega and appropriate zero points
	  if(i < 3) col[i] = 1e3*pow(10,0.4*(zp[i]-col[i])); //convert to mJy
	  if(i >= 3) col[i] = 1e3*pow(10,0.4*(zp[i-3]-col[i])); //convert to mJy
	}
      }
      catch(Table::NoSuchColumn){
	printf("Column %s does not exist in %s\n",column[i].c_str(),fitsfile.c_str());
	exit(1);
      }
    }
    
    observations.reserve(tablesize);
    double fluxes[3];
    bool accepted;

    for (unsigned int i=0;i<tablesize;i++){
      accepted = true;
      for (int j=0;j<3;j++){
	if((col[j][i] < flim[j]) or isnan(col[j][i])){
	  accepted = false;
	  break;
	}
	fluxes[j]=col[j][i];
      }
      if(accepted){
	observations.push_back(make_shared<obs>(fluxes,axes));
      }
    }
  }
  catch(FITS::NoSuchHDU){
    printf("HDU %i does not exist inf %s\n",hdunum,fitsfile.c_str());
    exit(1);
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

void obs_lib::get_all_colors(vector<double> &c1, vector<double> &c2){
  c1.resize(get_snum());
  c2.resize(get_snum());
  for(unsigned int i=0;i<observations.size();i++)
    observations[i]->get_colors(c1[i],c2[i]);
}
