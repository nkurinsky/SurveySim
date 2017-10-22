#include "fracs.h"
#include "lumfunct.h"

fracs::fracs(int numTypes, string fitsfile, int sedflags[4]){
  if (numTypes <= 1) {
    this->numTypes = 1;
    generate = false;
    printf("No AGN Detected\n");
  } else {
    this->numTypes = numTypes;
    generate = true;
  }

  load_sedmix(fitsfile);

  _t1  = 0;
  _t2  = 0;
  _zbt = 0;
  lf = NULL;
}

void fracs::load_sedmix(string fitsfile) {
  // Open FITS file
  unique_ptr<FITS> pfits;
  try {
    pfits.reset(new FITS(fitsfile));
  } catch (CCfits::FITS::CantOpen) {
    cout << "Cannot open " << fitsfile << endl;
    exit(1);
  }

  HDU& info = pfits->pHDU();
  info.readAllKeys();
  auto settings = info.keyWord();

  int nmodels;
  if (settings.count("NTYPES") == 1) {
    settings["NTYPES"]->value(nmodels);
  } else {
    cout << "Keyword NTYPES not specified; exiting" << endl;
    exit(1);
  }

  auto extensions = pfits->extension();

  // Read in types
  vector<string> types;
  string type;
  for (int i = 1; i <= nmodels; i++) {
    char buffer[10];
    sprintf(buffer, "SEDTYPE%i", i);
    string field(buffer);
    if (settings.count(field) > 0) {
      settings[field]->value(type);
      types.push_back(type);
    }
    else{
      cout << "Keyword " << field << " not specified; exiting" << endl;
      exit(1);
    }
  }

  unsigned int lnum = extensions.begin()->second->numCols() - 1;

  // Store data from SED_MIX
  for (auto type = types.begin(); type != types.end(); type++) {
    string t = toLower(*type);
    if(t == "agn")
      fracData[t] = 0.25; //default AGN fraction to 0.25
    else
      fracData[t] = 0.0; //otherwise allow for 0.0 default
    extensions.find("LUMDATA"/*"SED_MIX"*/)->second->column(*type).read(sedmix[t], 0, lnum);
  }
}
 
void fracs::set_lumfunct(lumfunct *lf) {
  if (lf != NULL)
    this->lf = lf;
  else
    cout << "ERROR: NULL Pointer Passed to Simulator" << endl;
}

//set fractions and power law parameters
void fracs::set_params(double lpars[]) {
  set_frac (lpars[LF::parameter::fa0],"agn");
  set_frac (lpars[LF::paramter::fcom],"com");
  set_t1  (lpars[LF::parameter::t1]);
  set_t2  (lpars[LF::parameter::t2]);
  set_zbt (lpars[LF::parameter::zbt]);
}

void fracs::set_t1(double t1) {
  if (t1 != _t1) {
    _t1=t1;
    evolZ.clear();
  }
}

void fracs::set_t2(double t2) {
  if (t2 != _t2) {
    _t2 = t2;
    evolZ.clear();
  }
}

void fracs::set_zbt(double zbt) {
  if (zbt != _zbt) {
    _zbt = zbt;
    evolZ.clear();
  }
}

void fracs::set_frac(double frac, string type) {
  if (frac < 0.0)
    frac = 0.0;
  else if (frac > 1.0)
    frac = 1.0;

  fracData[toLower(type)] = frac;
}

void fracs::set_fracs(unordered_map<string,double> newFracs) {
  // Maintain old info that isn't set by newFracs
  for (auto itr : newFracs) {
    fracData[toLower(itr.first)] = itr.second;
  }
}

//these return the fractions based on
// fractions for type,L and z read from fitsfile (which is at z=0)
// x redshift evolution based on fa0, t1, t2, and zbt parameters
// and fraction determined in fracData array
double fracs::get_frac(double lum, double redshift, string type) {
  // here fagn0 gives the fraction of AGN at L_ir=10^12 at z=0
  // zbt gives the break redshift for the sed type
  // t1 and t2 give the power relations below and after the break respectively

  if (evolZ.count(redshift) == 0) {
    if (redshift <= _zbt)
      evolZ[redshift] = pow((1 + redshift), _t1);
    if (redshift > _zbt)
      evolZ[redshift] = pow((1 + (redshift-_zbt)), _t1) * pow((1 + redshift), _t2);
  }

  //the z=0 fraction as read from the fits file
  //default assumption at L=10.^12 is of fSFG=0.75,fAGN0=0.25, fCOM=0
  double l = sedmix[type][(int)round(lum)];
  //Caution! for now only call this for type AGN
  //since that's the one that is subject to this evolZ
  double frac = (fracData[toLower(type)]/0.25) * l * evolZ[redshift];
  
  return (frac > 1 ? 1.0 : frac);
}

int fracs::get_sedtype(double lum, double redshift,int sedflags[NSEDS]) {
  // This section needs to be changed if any other mix
  // currently the SED templates are: SFG, COM, AGN, and COLD in that order
  double fracArray[numTypes];

  //Chase's code -- but added explicit suppression of templates if not used
  float fagn   = get_frac(lum, redshift, "agn");
  if(sedflags[2] == 0) fagn = 0;
  fracArray[0] = 1 - fagn; //SFG
  fracArray[1] = fracData["com"]  * fagn; //COMPOSITES are a fraction of the AGN
  if(sedflags[1] == 0) fracArray[1] = 0;  
  fracArray[2] = fagn;     //AGN
  
  if (!generate) return 0;

  float rnd = rng.flat(0, 1);
  for (int i = 0; i < numTypes; i++) {
    if (rnd < fracArray[i]) return i;
    else rnd -= fracArray[i];
  }
  
  return 0;
}
