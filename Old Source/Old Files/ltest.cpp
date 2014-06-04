#include "simulator.h"
#include "functions.h"
#include "cosmo.h"

using namespace std;

int main(){
  string modfile("model.fits");

  FITS *pInfile;
  pInfile = new FITS(modfile,Read);
  
  //reading primary header from fits file
  HDU& params_table = pInfile->pHDU();

//=================================================================  
//Read-in Luminosity Function Parameters
//-----------------------------------------------------------------

//this array holds phi0,L0,alpha,beta,p, and q
  double lpars[6];
  params_table.readKey("PHI0",lpars[0]);
  params_table.readKey("L0",lpars[1]);
  params_table.readKey("ALPHA",lpars[2]);
  params_table.readKey("BETA",lpars[3]);
  params_table.readKey("P",lpars[4]);
  params_table.readKey("Q",lpars[5]);
  
  lpars[4] = 6.0;
  lpars[1] = 10.0;

  delete pInfile;

  lumfunct lf;
  lf.set_phi0(lpars[0]);
  lf.set_L0(lpars[1]);
  lf.set_alpha(lpars[2]);
  lf.set_beta(lpars[3]);
  lf.set_p(lpars[4]);
  lf.set_q(lpars[5]);


  for (int i=0;i<6;i++)
    printf("\nLpar%i: %f",i,lpars[i]);

  double area, dz;
  int nz;
  area=pow((1.4*M_PI/180.0),2);
  nz = 10;
  dz = 0.5;

  double lums[] = {9,9.25,9.5,9.75,10.0,10.25,10.5,10.75,11.0,11.25,11.5,11.75,12.0,12.25,12.5,12.75,13.0};

  double zarray[nz];
  for (int i=0;i<nz;i++)
    zarray[i] = 0.1+i*dz;

  unsigned long nsrcs[nz][17];
  double tmpz,vol,high,scale[nz];
  for (int i=0;i<nz;i++){
    tmpz = zarray[i]+dz/2.0;
    vol = (dvdz(tmpz,area)*dz);
    high = lf.get_nsrcs(zarray[i],10.0)*vol;
    scale[i] = 10000e0/high;
    cout << "\nz = " << zarray[i] << ", s = " << scale[i] << ", high = " << high << endl;
    cout << "N: ";
    for(int j=0;j<17;j++){
      nsrcs[i][j] = (unsigned long)(scale[i]*vol*lf.get_nsrcs(zarray[i],lums[j]));
      cout << "{" << lums[j] << ", " << nsrcs[i][j] << "}, ";
    }
  }
  return 0;
}
