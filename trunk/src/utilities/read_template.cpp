//standard includes used by all programs in simulation
#include <cstdio>
#include <stdio.h>
#include <string>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <memory>
#include <cctype>

//For FITS input/output
#include <CCFits/CCfits>

using namespace CCfits;
using namespace std;

int main(){

  string fitsfile("Kirkpatrick2015_templates.fits");

  unique_ptr<FITS> pfits;
  try{
    pfits.reset(new FITS(fitsfile));
  }
  catch(CCfits::FITS::CantOpen){
    cout << "Cannot open "<< fitsfile << endl;
    return 1;
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
  double zmin;
  for(int i=0;i<nzbins;i++){
    char buffer[10];
    sprintf(buffer,"ZMIN%i",i);
    string field(buffer);
    if(settings.count(field) > 0){
      settings[field]->value(zmin);
    }
    else{
      printf("Keyword %s missing from header (%s)\n",field.c_str(),fitsfile.c_str());
      exit(1);
    }
    zmins.push_back(zmin);    
  }

  //Read redshift bins
  double zmax;
  for(int i=0;i<nzbins;i++){
    char buffer[10];
    sprintf(buffer,"ZMAX%i",i);
    string field(buffer);
    if(settings.count(field) > 0){
      settings[field]->value(zmax);
    }
    else{
      printf("Keyword %s missing from header (%s)\n",field.c_str(),fitsfile.c_str());
      exit(1);
    }
    zmaxs.push_back(zmax);
  }

  auto extensions=pfits->extension();

  //read luminosities from primary header
  int lnum=extensions.begin()->second->numCols()-1;
  unique_ptr<double[]> lums(new double[lnum]);
  unique_ptr<int[]> inds(new int[lnum]);
  char buffer[10];
  double temp;
  for (unsigned int i=0;i<lnum;i++) {
    std::string a = "L";
    sprintf(buffer,"%i",i+1);
    a.append(buffer);
    if(settings.count(a) == 1){
      settings[a]->value(temp);
      lums[i] = temp;
      inds[i] = i;
    }
    else{
      printf("Keyword %s missing from header (%s)\n",a.c_str(),fitsfile.c_str());
      exit(1);
    }
  }

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

  for (auto itr=type_exts.begin(); itr != type_exts.end(); ++itr){
    int type_len=itr->first.length();
    vector<double> lambda;
    vector<double> zs;
    for(auto eitr=itr->second.begin();eitr!=itr->second.end();++eitr){
      string extstr(*eitr);
      extstr.erase(0,type_len+2);
      int zbin(atoi(extstr.c_str()));
      zs.push_back(zmins[zbin]);
    }

    auto firstInstance=extensions.find(itr->second.front())->second;
    int tablelength=firstInstance->rows();
    firstInstance->column(1).read(lambda,0,tablelength);
    static double scale_factor(pow(10,scale));
    for(auto lam=lambda.begin();lam!=lambda.end();++lam)
      *lam *= scale_factor;
    cout << itr->first << "\t" << lambda.front() << " -> " << lambda.back() << endl;

    for(int i=2;i<lnum+2;i++){
      vector<double> fluxtemp;
      vector<double> fluxes;
      for(auto eitr=itr->second.begin();eitr!=itr->second.end();eitr++){
	extensions.find(*eitr)->second->column(i).read(fluxtemp,0,tablelength);
	for(int i=0;i<fluxtemp.size();i++){
	  fluxes.push_back(fluxtemp[i]);
	}
      }
      //initialize SED
    }
  }

  return 0;
}
