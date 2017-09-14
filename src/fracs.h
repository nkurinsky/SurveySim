// -*-c++-*-

#ifndef FRACS_H
#define FRACS_H

#include <math.h>
#include "lumfunct.h"
#include "functions.h"

using namespace CCfits;

class fracs {
  private:
    lumfunct *lf;
    int  numTypes;
    bool generate;
    RandomNumberGenerator rng;
    map <tuple<double,double>,double> values;
    map <double,double> lumPower;
    map <double,double> evolZ;
    double _t1;
    double _t2;
    double _zbt;
    unordered_map<string, double>         fracData;
    unordered_map<string, vector<double>> sedmix;
  public:
    fracs(int types, string fitsfile);
    // Setters
    void set_lumfunct(lumfunct *lf);
    void set_params(double lpars[]);
    void set_t1(double t1);
    void set_t2(double t2);
    void set_zbt(double zbt);
    void set_frac(double frac, string type);
    void set_fracs(unordered_map<string,double> fracs);
    void load_sedmix(string fitsfile);
    // Getters                // Defined here
    double get_t1()              {return _t1;}
    double get_t2()              {return _t2;}
    double get_zbt()             {return _zbt;}
    double get_frac(string type) {return fracData[type];}
    double get_frac(double lum, double redshift, string type);
    int get_sedtype(double lum, double redshift);
};


#endif
