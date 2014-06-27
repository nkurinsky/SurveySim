#include "filters.h"

using namespace std;

int main(){
  
  filter_lib my_filters("/Users/noahkurinsky/SurveySim/Source/filters/filterlib.txt");
  my_filters.load_filter(0,"SPIRE_250");
  my_filters.load_filter(1,"SPIRE_350");
  my_filters.load_filter(2,"SPIRE_500");

  my_filters.get(2).print(false);
  my_filters.get(0).print(false);
  my_filters.get(1).print(true);

  printf("\n\t%lf\n",my_filters.get(0).transmission(250e-6));
  printf("\t%lf\n",my_filters.get(1).transmission(250e-6));
  printf("\t%lf\n",my_filters.get(2).transmission(250e-6));

  return 0;
}
