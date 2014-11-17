AUTOMAKE_OPTION=foreign
SUBDIRS = alg_lib

AM_CXXFLAGS=-Wall -I/usr/local/include -I./ -std=c++11 -L ./alg_lib/

bin_PROGRAMS = fitter
fitter_SOURCES =  constants.h cosmo.h functions.h filters.h hist_lib.h lumfunct.h mc_util.h obs_lib.h sed_lib.h simulator.h simulator_utils.h numberCounts.h main.cpp cosmo.cpp functions.cpp filters.cpp hist_lib.cpp lumfunct.cpp mc_util.cpp obs_lib.cpp sed_lib.cpp simulator.cpp simulator_utils.cpp numberCounts.cpp
fitter_LDADD = libalglib.a
