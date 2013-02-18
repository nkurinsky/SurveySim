g++ -c fit_mcmc.cpp simulator.cpp lumfunct.cpp functions.cpp hist_lib.cpp model_lib.cpp 
#g++ test_gauss.cpp -g -lgsl -o test_gauss
g++ -g -lgsl -/scisoft/CCfits/include/CCfits -/scisoft/cfitsio/lcfitsio -lccfits fit_mcmc.o simulator.o lumfunct.o hist_lib.o model_lib.o -o fitter #-Wall -Werror -fexceptions -lgslcblas -O3
