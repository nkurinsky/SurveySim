AUTOMAKE_OPTION=foreign

CFLAGS=-Wall -I/usr/local/include -I./
AM_LDFLAGS=-L./ -lm -lgsl -lgslcblas -lcfitsio -lCCFits

bin_PROGRAMS = fitter
fitter_SOURCES =  constants.h cosmo.h functions.h hist_lib.h lumfunct.h mc_util.h obs_lib.h sed_lib.h simulator.h main.cpp functions.cpp hist_lib.cpp lumfunct.cpp mc_util.cpp obs_lib.cpp sed_lib.cpp simulator.cpp
