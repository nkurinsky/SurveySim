#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])

AC_INIT([SurveySim], [1.0], [nkurinsky91@gmail.com])
AC_CONFIG_SRCDIR([src/constants.h])

AM_MAINTAINER_MODE([disable])
AM_SILENT_RULES([yes])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_RANLIB

# Checks for libraries.

AC_CHECK_LIB([m],[cos],[],[
	echo "Error: LM library required, not found"; exit 1])
AC_CHECK_LIB([gslcblas],[cblas_sgemm],[],[
        echo "Error: GSLCBLAS library required, not found"; exit 1])
AC_CHECK_LIB([gsl],main,[],[
	echo "Error: GSL library required, not found"; exit 1])
AC_CHECK_LIB([cfitsio],[ffopen],[],[
	echo "Error:cfitsio library required, not found"; exit 1])
AC_CHECK_LIB([CCfits],[main],[],[
	echo "Error:CCFits library required, not found"; exit 1])

# Checks for header files.
AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([stdlib.h])
AC_CHECK_HEADERS([fitsio.h],[],[
	echo "Error: fitsio.h header file required, not found"; exit 1])
AC_CHECK_HEADERS([CCfits/CCfits.h],[],[
        echo "Error: CCfits/CCfits.h header file required, not found"; exit 1])
AC_CHECK_HEADERS([gsl/gsl_rng.h],[],[
	echo "Error: gls/gsl_rng.h header file required, not found"; exit 1])
AC_LANG_POP([C++])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt])

AC_CONFIG_FILES([Makefile src/Makefile filters/Makefile templates/Makefile lib_aux/Makefile src/alg_lib/Makefile python/Makefile python/SurveySim/Makefile examples/Makefile])
AC_OUTPUT
