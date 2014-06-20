#include "filters.h"

using namespace std;

int main(){
  
  filter_lib my_filters("/Users/noahkurinsky/Desktop/filt/filterlist.txt");
  my_filters.load_filter(0,"SPIRE_250");
  my_filters.load_filter(1,"SPIRE_350");
  my_filters.load_filter(2,"SPIRE_500");

  my_filters.get(1).print(true);

  printf("\n%lf\n",my_filters.get(1).transmission(3700000.0));

  return 0;
}
