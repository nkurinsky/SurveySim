#!/bin/bash

printf "Building Simulator.cpp..."
g++ -c simulator.cpp 
if [ "$?" = "0" ]; then
    printf "Successful \n"
else
    printf "\nSimulator.cpp build failed...exiting\n\n" 1>&2
    exit 1
fi

printf "Building Lumfunct.cpp..." 
g++ -c lumfunct.cpp 

g++ -c functions.cpp 
g++ -c hist_lib.cpp 
g++ -c sed_lib.cpp 

echo "
main build
"

g++ -g -lgsl -lcfitsio -lccfits fit_mcmc.cpp simulator.o lumfunct.o hist_lib.o sed_lib.o -o fitter #-Wall -Werror -fexceptions -lgslcblas -O3
