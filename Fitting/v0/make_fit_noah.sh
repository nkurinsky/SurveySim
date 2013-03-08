#!/bin/bash

clear
printf "\nSurvey Simulator Compilation Code\n"
printf "Last Updated 4/7/13 by Noah Kurinsky\n\n"

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
if [ "$?" = "0" ]; then
    printf "Successful \n"
else
    printf "\nLumfunct.cpp build failed...exiting\n\n" 1>&2
    exit 1
fi

printf "Building Functions.cpp..." 
g++ -c functions.cpp 
if [ "$?" = "0" ]; then
    printf "Successful \n"
else
    printf "\nFunctions.cpp build failed...exiting\n\n" 1>&2
    exit 1
fi

printf "Building hist_lib.cpp..." 
g++ -c hist_lib.cpp 
if [ "$?" = "0" ]; then
    printf "Successful \n"
else
    printf "\nhist_lib.cpp build failed...exiting\n\n" 1>&2
    exit 1
fi

printf "Building sed_lib.cpp..." 
g++ -c sed_lib.cpp 
if [ "$?" = "0" ]; then
    printf "Successful \n"
else
    printf "\nsed_lib.cpp build failed...exiting\n\n" 1>&2
    exit 1
fi

printf "Building obs_lib.cpp..."
g++ -c obs_lib.cpp
if [ "$?" = "0" ]; then
    printf "Successful \n"
else
    printf "\nobs_lib.cpp build failed...exiting \n\n" 1>&2
    exit 1
fi

printf "\nBuilding Main..."
g++ -g -lgsl -lcfitsio -lccfits fit_mcmc.cpp simulator.o lumfunct.o hist_lib.o sed_lib.o obs_lib.o functions.o -o fitter -Wall -Werror -fexceptions -lgslcblas -O3
if [ "$?" = "0" ]; then
    printf "Successful \n"
    exit 0
else
    printf "\nMain build failed...exiting\n\n" 1>&2
    exit 1
fi