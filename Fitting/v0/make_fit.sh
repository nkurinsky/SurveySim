g++ -c fit_mcmc.cpp simulator.cpp lumfunct.cpp hist_lib.cpp model_lib.cpp 
# g++ -g -lgsl -/scisoft/CCfits/include/CCfits -/scisoft/cfitsio/lcfitsio -lccfits simulator.o lumfunct.o hist_lib.o model_lib.o -o fitter -Wall -Werror -fexceptions -lgslcblas -O3
