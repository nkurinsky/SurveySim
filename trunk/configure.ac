#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])

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
LIBS="-lm -lgslcblas -lgsl -lcfitsio -lCCFits"

AC_CHECK_LIB(m,cos)
AC_CHECK_LIB([gslcblas],[cblas_sgemm],[],[
        echo "Warning: GSLCBLAS library required, not found"])
AC_CHECK_LIB([gsl],main,[],[
	echo "Warning: GSL library required, not found"])
AC_CHECK_LIB([cfitsio],[ffopen],[],[
	echo "Warning:cfitsio library required, not found"])
AC_CHECK_LIB([CCFits],[main],[],[
	echo "Warning:CCFits library required, not found"])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt])

AC_CONFIG_FILES([Makefile src/Makefile filters/Makefile templates/Makefile obs/Makefile lib_aux/Makefile src/alg_lib/Makefile])
AC_OUTPUT
