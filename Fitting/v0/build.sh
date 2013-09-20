#!/bin/bash
##################################################################
#
# "Fitter" compilation script, written by Noah Kurinsky
# Script detects host and changes main build settings accordingly
#
##################################################################
clear

printf "\nSurvey Simulator Compilation Code\n"
printf "Last Updated 9/20/13 by Noah Kurinsky\n\n"
printf "Building using settings for "$HOST"\n"

printf "Building Simulator.cpp..."
g++ -g -c simulator.cpp -Wall -Wextra -Werror
if [ "$?" = "0" ]; then
    printf "Successful \n"
else
    printf "\nSimulator.cpp build failed...exiting\n\n" 1>&2
    exit 1
fi

printf "Building Lumfunct.cpp..." 
g++ -g -c lumfunct.cpp -Wall -Werror -Wextra
if [ "$?" = "0" ]; then
    printf "Successful \n"
else
    printf "\nLumfunct.cpp build failed...exiting\n\n" 1>&2
    exit 1
fi

printf "Building Functions.cpp..." 
g++ -g -c functions.cpp -Wall -Werror -Wextra
if [ "$?" = "0" ]; then
    printf "Successful \n"
else
    printf "\nFunctions.cpp build failed...exiting\n\n" 1>&2
    exit 1
fi

printf "Building hist_lib.cpp..." 
g++ -g -c hist_lib.cpp -Wall -Werror -Wextra
if [ "$?" = "0" ]; then
    printf "Successful \n"
else
    printf "\nhist_lib.cpp build failed...exiting\n\n" 1>&2
    exit 1
fi

printf "Building sed_lib.cpp..." 
g++ -g -c sed_lib.cpp -Wall -Werror -Wextra
if [ "$?" = "0" ]; then
    printf "Successful \n"
else
    printf "\nsed_lib.cpp build failed...exiting\n\n" 1>&2
    exit 1
fi

printf "Building obs_lib.cpp..."
g++ -g -c obs_lib.cpp -Wall -Werror -Wextra
if [ "$?" = "0" ]; then
    printf "Successful \n"
else
    printf "\nobs_lib.cpp build failed...exiting \n\n" 1>&2
    exit 1
fi

printf "\nBuilding Main (fit_mcmc.cpp)..."
if [ $HOST = "ariel" ]; then
    printf "\nUsing settings for Ariel..."
    g++ -g -lgsl -lcfitsio -lccfits fit_mcmc.cpp simulator.o lumfunct.o hist_lib.o sed_lib.o obs_lib.o functions.o -o fitter -Wall -Werror -Wextra -fexceptions -lgslcblas -O3
    if [ "$?" = "0" ]; then
	printf "Successful \n"
    else
	printf "\nMain build failed...exiting\n\n" 1>&2
	exit 1
    fi
else
    printf "\nUsing settings for ASajina..."
    g++ -g -lgsl -/scisoft/CCfits/include/CCfits -/scisoft/cfitsio/lcfitsio -lccfits fit_mcmc.cpp simulator.o lumfunct.o hist_lib.o sed_lib.o obs_lib.o functions.o -o fitter -Wall -Werror -Wextra -fexceptions -lgslcblas -O3
    if [ "$?" = "0" ]; then
	printf "Successful \n"
    else
	printf "\nMain build failed...exiting\n\n" 1>&2
    fi
fi

printf "\nProgram successfully compiled as \"fitter\"\n"

printf "Copying to Widget directory..."
cp fitter ../../Widget/
if [ "$?" = "0" ]; then
    printf "Successful \n"
else
    printf "failed for unknown reason, check directory structure"
fi

printf "\nBuild Successful...Exiting\n\n"

exit 0