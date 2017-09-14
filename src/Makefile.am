AUTOMAKE_OPTION=foreign
SUBDIRS = alg_lib

AM_CXXFLAGS=-Wall -I/usr/local/include -I./ -I./alg_lib -std=c++11 -g

bin_PROGRAMS = SurveySim
SurveySim_SOURCES =  SurveySim.h fracs.h constants.h cosmo.h functions.h filters.h hist_lib.h lumfunct.h mc_util.h obs_lib.h sed_lib.h simulator.h simulator_utils.h numberCounts.h SurveySim.cpp fracs.cpp main.cpp cosmo.cpp functions.cpp filters.cpp hist_lib.cpp lumfunct.cpp mc_util.cpp obs_lib.cpp sed_lib.cpp simulator.cpp simulator_utils.cpp numberCounts.cpp
SurveySim_LDADD = ./alg_lib/libalglib.a
