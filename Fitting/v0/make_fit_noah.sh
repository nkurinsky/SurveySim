g++ -c fit_mcmc.cpp simulator.cpp lumfunct.cpp functions.cpp hist_lib.cpp model_lib.cpp 

g++ -g -lgsl -lcfitsio -lccfits fit_mcmc.o simulator.o lumfunct.o hist_lib.o model_lib.o -o fitter #-Wall -Werror -fexceptions -lgslcblas -O3
