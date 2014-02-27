#!/bin/bash
##################################################################
#
# "Fitter" compilation script, written by Noah Kurinsky
# Script detects host and changes main build settings accordingly
#
##################################################################
clear

printf "\nSurvey Simulator Compilation Code\n"
printf "Last Updated 11/8/13 by Noah Kurinsky\n\n"
printf "Building using settings for "$HOST"\n"

files=( "Simulator.cpp" "Lumfunct.cpp" "Functions.cpp" "hist_lib.cpp" "sed_lib.cpp" "obs_lib.cpp" "mc_util.cpp" )

for file in "${files[@]}"
do
    printf "Building "$file"..."
    g++ -g -c $file -Wall -Wextra -Werror
    if [ "$?" = "0" ]; then
	printf "Successful \n"
    else
	printf "\n"$file" build failed...exiting\n\n" 1>&2
	exit 1
    fi
done

mainfile="main.cpp"
printf "\nBuilding Main ("$mainfile")..."
if [ $HOST = "ariel" ]; then
    printf "\nUsing settings for Ariel..."
    g++ -g -lgsl -lcfitsio -lccfits $mainfile ./*.o -o fitter -Wall -Werror -Wextra -fexceptions -lgslcblas -O3
    if [ "$?" = "0" ]; then
	printf "Successful \n"
    else
	printf "\nMain build failed...exiting\n\n" 1>&2
	exit 1
    fi
else
    printf "\nUsing settings for ASajina..."
    g++ -g -lgsl -/scisoft/CCfits/include/CCfits -/scisoft/cfitsio/lcfitsio -lccfits $mainfile *.o -o fitter -Wall -Werror -Wextra -fexceptions -lgslcblas -O3
    if [ "$?" = "0" ]; then
	printf "Successful \n"
    else
	printf "\nMain build failed...exiting\n\n" 1>&2
    fi
fi

printf "\nProgram successfully compiled as \"fitter\"\n"

printf "Copying to Widget directory..."
cp -r fitter* ../../Widget/
if [ "$?" = "0" ]; then
    printf "Successful \n"
else
    printf "failed for unknown reason, check directory structure"
fi

printf "\nBuild Successful\n\n"

rm *.o

exit 0