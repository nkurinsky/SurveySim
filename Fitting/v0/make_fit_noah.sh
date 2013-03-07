g++ -c simulator.cpp 
g++ -c lumfunct.cpp 
g++ -c functions.cpp 
g++ -c hist_lib.cpp 
g++ -c sed_lib.cpp 

echo "
main build
"

g++ -g -lgsl -lcfitsio -lccfits fit_mcmc.cpp simulator.o lumfunct.o hist_lib.o sed_lib.o -o fitter #-Wall -Werror -fexceptions -lgslcblas -O3
