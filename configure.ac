## Process this file with autoconf to produce a configure script.
##
## Configure.ac for Rlibeemd
##
## Modified from RcppGSL's (Francois & Eddelbuettel) example package

# The version set here will propagate to other files from here
AC_INIT([Rlibeemd],[1.4.2])

## Use gsl-config to find arguments for compiler and linker flags
##
## Check for non-standard programs: gsl-config(1)
AC_PATH_PROG([GSL_CONFIG], [gsl-config])

## If gsl-config was found, use it
if test "${GSL_CONFIG}" != ""; then
    # Use gsl-config for header and linker arguments (without BLAS which we get from R)
    GSL_CFLAGS=`${GSL_CONFIG} --cflags`
    GSL_LIBS=`${GSL_CONFIG} --libs`
else
    AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi

# Substitute these variables in src/Makevars.in used to create src/Makevars
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
