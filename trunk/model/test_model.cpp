#include <iostream>
#include <CCFits/CCfits>

using namespace CCfits;
using namespace std;

int main(){

  FITS *pInfile;
  string testfile("model.fits");
  pInfile=new FITS(testfile,Read,true);
  ExtHDU& filters = pInfile->extension(1);
  int ntcols,ntrows;
  ntcols=filters.numCols();
  if(ntcols != 6){
     printf("Wrong number of filters included (need 3):");
     exit(1);
   }
  //cout<<filters.column(1);
  return 0;
  }
